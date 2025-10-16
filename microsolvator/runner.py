"""High-level API to orchestrate CREST microsolvation runs."""

from __future__ import annotations

import os
import shlex
import subprocess
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Callable, Optional, Sequence, Set

from ase import Atoms
from ase.io import read as ase_read, write as ase_write

from .command import build_crest_command
from .config import MicrosolvatorConfig
from .install import resolve_crest_binary, resolve_xtb_binary
from .results import MicrosolvationResult
from .support import validate_implicit_choice


RunCommand = Callable[..., subprocess.CompletedProcess[str]]


class Microsolvator:
    """Coordinate CREST microsolvation jobs from ASE inputs."""

    @classmethod
    def run(
        cls,
        *,
        solute: Atoms,
        solvent: Atoms,
        config: MicrosolvatorConfig,
        constrained_indices: Optional[Sequence[int]] = None,
        constrain_solute: bool = False,
        keep_temps: bool = False,
        working_directory: Optional[Path] = None,
        run_command: Optional[RunCommand] = None,
        prepare_only: bool = False,
    ) -> MicrosolvationResult:
        """Execute a CREST microsolvation run and return parsed results."""

        if config.implicit_model and config.implicit_solvent:
            validate_implicit_choice(
                method=config.method,
                model=config.implicit_model,
                solvent=config.implicit_solvent,
            )

        executor = run_command or _default_runner

        if prepare_only and not keep_temps and working_directory is None:
            keep_temps = True

        if working_directory is not None:
            workdir = Path(working_directory)
            workdir.mkdir(parents=True, exist_ok=True)
            return cls._execute(
                solute=solute,
                solvent=solvent,
                config=config,
                workdir=workdir,
                constrained_indices=constrained_indices,
                constrain_solute=constrain_solute,
                executor=executor,
                prepare_only=prepare_only,
            )

        if keep_temps:
            workdir = Path(tempfile.mkdtemp(prefix="microsolvator_"))
            return cls._execute(
                solute=solute,
                solvent=solvent,
                config=config,
                workdir=workdir,
                constrained_indices=constrained_indices,
                constrain_solute=constrain_solute,
                executor=executor,
                prepare_only=prepare_only,
            )

        with tempfile.TemporaryDirectory(prefix="microsolvator_") as tmpdir:
            workdir = Path(tmpdir)
            return cls._execute(
                solute=solute,
                solvent=solvent,
                config=config,
                workdir=workdir,
                constrained_indices=constrained_indices,
                constrain_solute=constrain_solute,
                executor=executor,
                prepare_only=prepare_only,
            )

    @classmethod
    def _execute(
        cls,
        *,
        solute: Atoms,
        solvent: Atoms,
        config: MicrosolvatorConfig,
        workdir: Path,
        constrained_indices: Optional[Sequence[int]],
        constrain_solute: bool,
        executor: RunCommand,
        prepare_only: bool,
    ) -> MicrosolvationResult:
        workdir.mkdir(parents=True, exist_ok=True)

        solute_path = workdir / "solute.xyz"
        solvent_path = workdir / "solvent.xyz"

        ase_write(solute_path, solute, format='xyz')
        ase_write(solvent_path, solvent, format='xyz')

        constraints_written = _write_constraints(
            workdir=workdir,
            solute_length=len(solute),
            constrained_indices=constrained_indices,
            constrain_solute=constrain_solute,
        )

        exec_config = config
        if constraints_written and not config.nopreopt:
            exec_config = replace(config, nopreopt=True)

        crest_executable = resolve_crest_binary(exec_config.crest_executable)
        xtb_executable = resolve_xtb_binary(exec_config.xtb_executable)

        command = build_crest_command(
            config=exec_config,
            crest_executable=crest_executable,
            xtb_executable=xtb_executable,
            solute_path=solute_path,
            solvent_path=solvent_path,
        )

        if prepare_only:
            command_str = shlex.join(command)
            print(command_str)
            return MicrosolvationResult(
                command=command,
                working_directory=workdir,
                best_structure=None,
                ensemble=[],
                population_path=None,
                stdout="",
                stderr="",
                executed=False,
            )

        env = _build_subprocess_env(
            xtb_executable=command[command.index("--xnam") + 1],
            crest_executable=command[0],
        )

        try:
            completed = executor(command, workdir, env)
        except TypeError:
            completed = executor(command, workdir)

        best_structure = _read_optional_atoms(workdir / "crest_best.xyz")
        ensemble = _read_optional_ensemble(workdir / "full_ensemble.xyz")
        population_path = (workdir / "full_population.dat") if (workdir / "full_population.dat").exists() else None

        result = MicrosolvationResult(
            command=command,
            working_directory=workdir,
            best_structure=best_structure,
            ensemble=ensemble,
            population_path=population_path,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
        result.ensure_outputs()
        return result


def _default_runner(
    command: Sequence[str],
    workdir: Path,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=str(workdir),
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )


def _build_subprocess_env(
    *,
    xtb_executable: str,
    crest_executable: str,
) -> dict[str, str]:
    env = os.environ.copy()
    env["CREST_BIN"] = crest_executable
    env["XTB_BIN"] = xtb_executable
    return env


def _write_constraints(
    *,
    workdir: Path,
    solute_length: int,
    constrained_indices: Optional[Sequence[int]],
    constrain_solute: bool,
) -> bool:
    indices: Set[int] = set()

    if constrain_solute:
        indices.update(range(1, solute_length + 1))

    if constrained_indices:
        for index in constrained_indices:
            if index <= 0:
                raise ValueError("Constraint indices must be 1-based and positive")
            indices.add(index)

    if not indices:
        return False

    content_lines = ["$constrain", "atoms: " + ",".join(str(i) for i in sorted(indices)), "$end", ""]
    (workdir / ".xcontrol").write_text("\n".join(content_lines), encoding="utf-8")
    return True


def _read_optional_atoms(path: Path) -> Optional[Atoms]:
    if not path.exists():
        return None
    try:
        return ase_read(path)
    except FileNotFoundError:
        return None


def _read_optional_ensemble(path: Path) -> list[Atoms]:
    if not path.exists():
        return []
    try:
        return list(ase_read(path, index=":"))
    except FileNotFoundError:
        return []
