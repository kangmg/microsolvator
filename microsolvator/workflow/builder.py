"""Main orchestrator for the solvated reaction trajectory workflow."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple

from ase import Atoms

from ..results import MicrosolvationResult
from ..runner import Microsolvator
from .config import (
    EquilibrationConfig,
    MicrosolvatorConfig,
    PackmolConfig,
    RelaxationConfig,
    SolvationWorkflowConfig,
)
from .equilibration import equilibrate
from .packmol import PackmolSolvator
from .results import SolvatedTrajectoryResult
from .swap import relax_interface, swap_solute


class SolvatedTrajectoryBuilder:
    """Build explicitly solvated reaction trajectories.

    The full pipeline is::

        microsolvate_ts  →  pack_solvent_box  →  equilibrate  →  swap_and_relax

    Each step is also exposed as an independent ``@staticmethod`` so that
    intermediate results can be reused or inspected without re-running the
    whole workflow.

    Example
    -------
    >>> result = SolvatedTrajectoryBuilder.build(
    ...     reaction_images=images,
    ...     solvent=water,
    ...     config=config,
    ...     calculator=calc,
    ... )
    >>> write("solvated_guess.traj", result.solvated_images)
    """

    @classmethod
    def build(
        cls,
        *,
        reaction_images: List[Atoms],
        solvent: Atoms,
        config: SolvationWorkflowConfig,
        calculator,
        log_callback: Optional[Callable[[str, Any], None]] = None,
    ) -> SolvatedTrajectoryResult:
        """Run the full solvation workflow.

        Parameters
        ----------
        reaction_images:
            Initial NEB / string-method trajectory.  Must contain at least
            one structure.
        solvent:
            One solvent molecule.
        config:
            Full workflow configuration.
        calculator:
            ASE-compatible calculator used for MD and relaxation (Steps 3–4).
            Charge and multiplicity should be pre-configured on this object.
        log_callback:
            Optional hook called at each step boundary with
            ``(step_name, info_dict)``.

        Returns
        -------
        SolvatedTrajectoryResult
        """
        if not reaction_images:
            raise ValueError("reaction_images must contain at least one structure")

        ts_index = config.ts_index
        if ts_index is None:
            ts_index = len(reaction_images) // 2
        if not (0 <= ts_index < len(reaction_images)):
            raise ValueError(
                f"ts_index {ts_index} out of range "
                f"for {len(reaction_images)} images"
            )

        workdir = _resolve_workdir(config)

        # ── Step 1: Microsolvation ──────────────────────────────────────────
        if log_callback:
            log_callback("microsolvation", {"ts_index": ts_index})

        ts_image = reaction_images[ts_index]
        microsolv_result = cls.microsolvate_ts(
            ts_image=ts_image,
            solvent=solvent,
            config=config.microsolv,
            working_directory=workdir / "microsolvation" if workdir else None,
        )

        # ensemble=True  → best_structure (more refined)
        # ensemble=False → final (grow result only)
        microsolvated_ts = microsolv_result.best_structure or microsolv_result.final
        if microsolvated_ts is None:
            raise ValueError(
                "Microsolvation produced no valid cluster structure. "
                "Check CREST output in the working directory."
            )

        n_solute = len(ts_image)
        n_cluster = len(microsolvated_ts)

        # ── Step 2: Packmol box solvation ───────────────────────────────────
        if log_callback:
            log_callback("packmol", {"n_cluster": n_cluster})

        boxed_system, _ = cls.pack_solvent_box(
            cluster=microsolvated_ts,
            solvent=solvent,
            config=config.packmol,
            workdir=workdir / "packmol" if workdir else None,
        )

        # Propagate charge / multiplicity from the TS image
        charge = ts_image.info.get("charge", config.microsolv.charge)
        multiplicity = ts_image.info.get("multiplicity", config.microsolv.uhf + 1)
        boxed_system.info["charge"] = charge
        boxed_system.info["multiplicity"] = multiplicity

        # ── Step 3: MD equilibration ────────────────────────────────────────
        # Only freeze the solute atoms; the microsolvation shell is free to
        # equilibrate together with the bulk solvent so that the
        # shell–bulk interface relaxes naturally.
        if log_callback:
            log_callback("equilibration", {"n_fixed": n_solute})

        traj_path = str(workdir / "equilibration.traj") if workdir else None
        equilibrated_system = cls.equilibrate(
            system=boxed_system,
            n_fixed=n_solute,
            calculator=calculator,
            config=config.equilibration,
            traj_path=traj_path,
        )

        # ── Step 4: Swap endpoints & relax ──────────────────────────────────
        if log_callback:
            log_callback("swap_relax", {"n_images": len(reaction_images)})

        solvated_images = cls.swap_and_relax(
            template=equilibrated_system,
            reaction_images=reaction_images,
            n_solute=n_solute,
            ts_index=ts_index,
            calculator=calculator,
            config=config.relaxation,
        )

        # Propagate cell / PBC from the equilibrated template to all images
        for img in solvated_images:
            img.cell = equilibrated_system.cell.copy()
            img.pbc = equilibrated_system.pbc.copy()
            img.info["charge"] = charge
            img.info["multiplicity"] = multiplicity

        return SolvatedTrajectoryResult(
            solvated_images=solvated_images,
            microsolvated_ts=microsolvated_ts,
            boxed_system=boxed_system,
            equilibrated_system=equilibrated_system,
            n_solute_atoms=n_solute,
            n_cluster_atoms=n_cluster,
            n_total_atoms=len(equilibrated_system),
            ts_index=ts_index,
            config=config,
            working_directory=workdir or Path("."),
        )

    # ── Individual step methods ─────────────────────────────────────────────

    @staticmethod
    def microsolvate_ts(
        ts_image: Atoms,
        solvent: Atoms,
        config: MicrosolvatorConfig,
        **kwargs,
    ) -> MicrosolvationResult:
        """Step 1: CREST QCG microsolvation of the TS image.

        The solute geometry is always constrained (``constrain_solute=True``).
        All extra keyword arguments are forwarded to :meth:`Microsolvator.run`.
        """
        return Microsolvator.run(
            solute=ts_image,
            solvent=solvent,
            config=config,
            constrain_solute=True,
            **kwargs,
        )

    @staticmethod
    def pack_solvent_box(
        cluster: Atoms,
        solvent: Atoms,
        config: PackmolConfig,
        workdir: Optional[Path] = None,
    ) -> Tuple[Atoms, int]:
        """Step 2: Packmol box solvation.

        Returns ``(boxed_system, n_bulk_solvent)``.
        """
        return PackmolSolvator.run(
            cluster=cluster,
            solvent=solvent,
            config=config,
            workdir=workdir,
        )

    @staticmethod
    def equilibrate(
        system: Atoms,
        n_fixed: int,
        calculator,
        config: EquilibrationConfig,
        traj_path: Optional[str] = None,
        log_callback: Optional[Callable[[str, Any], None]] = None,
    ) -> Atoms:
        """Step 3: Constrained MD equilibration."""
        return equilibrate(
            system=system,
            n_fixed=n_fixed,
            calculator=calculator,
            config=config,
            traj_path=traj_path,
            log_callback=log_callback,
        )

    @staticmethod
    def swap_and_relax(
        template: Atoms,
        reaction_images: List[Atoms],
        n_solute: int,
        ts_index: int,
        calculator,
        config: RelaxationConfig,
    ) -> List[Atoms]:
        """Step 4: Swap each image's solute and relax the solvent.

        Propagates outward from the TS image in both directions so that
        each image inherits the relaxed solvent from the previous
        (closer-to-TS) image, ensuring continuity along the path.

        Only the solute atoms (indices ``0:n_solute``) are frozen during
        relaxation.  The microsolvation shell and bulk solvent are free
        to adapt to each image's solute geometry.

        The TS image (``ts_index``) is returned as a copy of *template*
        without any swap or relaxation.
        """
        n_images = len(reaction_images)
        results: List[Optional[Atoms]] = [None] * n_images

        # TS image: direct copy of the equilibrated template
        results[ts_index] = template.copy()

        # Forward: TS+1, TS+2, … → product
        current = template
        for i in range(ts_index + 1, n_images):
            swapped = swap_solute(current, reaction_images[i], n_solute)
            relaxed = relax_interface(swapped, n_solute, calculator, config)
            results[i] = relaxed
            current = relaxed

        # Backward: TS-1, TS-2, … → reactant
        current = template
        for i in range(ts_index - 1, -1, -1):
            swapped = swap_solute(current, reaction_images[i], n_solute)
            relaxed = relax_interface(swapped, n_solute, calculator, config)
            results[i] = relaxed
            current = relaxed

        return results  # type: ignore[return-value]


def _resolve_workdir(config: SolvationWorkflowConfig) -> Optional[Path]:
    """Return a working directory Path, or None for fully in-memory execution."""
    if config.working_directory is not None:
        p = Path(config.working_directory)
        p.mkdir(parents=True, exist_ok=True)
        return p
    if config.keep_temps:
        p = Path(tempfile.mkdtemp(prefix="solvated_traj_"))
        return p
    return None
