"""Packmol wrapper for box solvation (Step 2)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from ase import Atoms
from ase.io import read as ase_read, write as ase_write

from .config import PackmolConfig
from .utils import estimate_box_size, count_solvent_molecules, validate_packmol_output


def resolve_packmol(executable: Optional[str]) -> str:
    """Return the path to the packmol binary.

    Raises
    ------
    FileNotFoundError
        If the binary cannot be located.
    """
    name = executable or "packmol"
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(
            f"packmol executable '{name}' not found in PATH. "
            "Install packmol or set PackmolConfig.packmol_executable."
        )
    return path


def _write_packmol_input(
    input_path: Path,
    output_path: Path,
    cluster_path: Path,
    solvent_path: Path,
    box_center: np.ndarray,
    box_size: float,
    n_solvent: int,
    min_distance: float,
) -> None:
    """Write a Packmol input file.

    The cluster is placed with its geometric center at *box_center* using
    Packmol's ``center fixed`` directive (no rotation applied).
    Solvent molecules are packed uniformly inside the full box, subject to
    the *min_distance* tolerance constraint.
    """
    cx, cy, cz = box_center
    gap = min_distance
    upper = box_size - gap

    lines = [
        f"tolerance {min_distance:.4f}",
        "filetype xyz",
        f"output {output_path.as_posix()}",
        "",
        # Cluster: fixed at box centre, no rotation
        f"structure {cluster_path.as_posix()}",
        "  number 1",
        "  center",
        f"  fixed {cx:.6f} {cy:.6f} {cz:.6f} 0. 0. 0.",
        "end structure",
        "",
        # Bulk solvent: packed uniformly inside the box
        f"structure {solvent_path.as_posix()}",
        f"  number {n_solvent}",
        f"  inside box {gap:.4f} {gap:.4f} {gap:.4f} "
        f"{upper:.4f} {upper:.4f} {upper:.4f}",
        "end structure",
        "",
    ]
    input_path.write_text("\n".join(lines), encoding="utf-8")


class PackmolSolvator:
    """Run Packmol to solvate a cluster inside a cubic box."""

    @classmethod
    def run(
        cls,
        *,
        cluster: Atoms,
        solvent: Atoms,
        config: PackmolConfig,
        workdir: Optional[Path] = None,
    ) -> Tuple[Atoms, int]:
        """Solvate *cluster* and return ``(boxed_system, n_bulk_solvent)``.

        The returned ``Atoms`` object has atom ordering
        ``[cluster | bulk_solvent_mol_1 | … | bulk_solvent_mol_N]`` and its
        ``cell`` / ``pbc`` attributes set to the cubic simulation box.

        Parameters
        ----------
        cluster:
            Microsolvated TS cluster (solute + first shell).
        solvent:
            One solvent molecule.
        config:
            Packmol configuration.
        workdir:
            Directory for Packmol files; None → system temp (cleaned up).

        Returns
        -------
        tuple[Atoms, int]
            The solvated system and the number of bulk solvent molecules added.
        """
        if config.solvent_density is None and config.n_bulk_solvent is None:
            raise ValueError(
                "PackmolConfig: provide either solvent_density or n_bulk_solvent"
            )
        if config.solvent_density is not None and config.n_bulk_solvent is not None:
            raise ValueError(
                "PackmolConfig: provide solvent_density or n_bulk_solvent, not both"
            )

        packmol_bin = resolve_packmol(config.packmol_executable)
        margin = config.box_margin

        for attempt in range(config.max_retries + 1):
            try:
                return cls._attempt(
                    cluster=cluster,
                    solvent=solvent,
                    config=config,
                    margin=margin,
                    packmol_bin=packmol_bin,
                    workdir=workdir,
                )
            except subprocess.CalledProcessError:
                if attempt == config.max_retries:
                    raise
                margin += config.retry_margin_increment

        # unreachable
        raise RuntimeError("Packmol failed after all retries")

    @classmethod
    def _attempt(
        cls,
        *,
        cluster: Atoms,
        solvent: Atoms,
        config: PackmolConfig,
        margin: float,
        packmol_bin: str,
        workdir: Optional[Path],
    ) -> Tuple[Atoms, int]:
        use_temp = workdir is None
        tmpdir_str = tempfile.mkdtemp(prefix="packmol_") if use_temp else str(workdir)
        wd = Path(tmpdir_str)
        wd.mkdir(parents=True, exist_ok=True)

        try:
            box_size = estimate_box_size(cluster, margin)
            box_center = np.full(3, box_size / 2.0)

            if config.n_bulk_solvent is not None:
                n_solvent = config.n_bulk_solvent
            else:
                n_solvent = count_solvent_molecules(
                    box_size, solvent, config.solvent_density,  # type: ignore[arg-type]
                    cluster=cluster,
                )

            cluster_path = wd / "cluster.xyz"
            solvent_path = wd / "solvent.xyz"
            output_path = wd / "solvated.xyz"
            input_path = wd / "packmol.inp"

            ase_write(cluster_path, cluster, format="xyz")
            ase_write(solvent_path, solvent, format="xyz")

            _write_packmol_input(
                input_path=input_path,
                output_path=output_path,
                cluster_path=cluster_path,
                solvent_path=solvent_path,
                box_center=box_center,
                box_size=box_size,
                n_solvent=n_solvent,
                min_distance=config.min_distance,
            )

            with open(input_path, encoding="utf-8") as fh:
                subprocess.run(
                    [packmol_bin],
                    stdin=fh,
                    check=True,
                    capture_output=True,
                    text=True,
                )

            boxed = ase_read(output_path)
            validate_packmol_output(boxed, cluster, n_solvent, len(solvent))

            boxed.cell = [box_size, box_size, box_size]
            boxed.pbc = True

            return boxed, n_solvent

        finally:
            if use_temp:
                shutil.rmtree(tmpdir_str, ignore_errors=True)
