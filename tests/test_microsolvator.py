from pathlib import Path

import subprocess

from ase import Atoms
from ase.io import write as ase_write

from microsolvator import (
    Microsolvator,
    MicrosolvatorConfig,
    list_supported_implicit_solvents,
    supports_implicit_solvent,
)


def test_support_tables_consistency():
    assert supports_implicit_solvent(method="gfn2", model="alpb", solvent="h2o")
    data = list_supported_implicit_solvents(method="gfn2", model="alpb")
    assert "gfn2" in data and "alpb" in data["gfn2"]
    assert "h2o" in data["gfn2"]["alpb"]


def test_microsolvator_run_with_constraints_and_mock_executor(tmp_path: Path):
    solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
    solvent = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])

    config = MicrosolvatorConfig(
        nsolv=3,
        method="gfn2",
        implicit_model="alpb",
        implicit_solvent="h2o",
    )

    def fake_runner(command, workdir):
        best = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.8, 0]])
        ensemble = [best, best.copy()]
        ase_write(workdir / "crest_best.xyz", best)
        ase_write(workdir / "full_ensemble.xyz", ensemble)
        (workdir / "full_population.dat").write_text("1 0.5\n", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="crest ok",
            stderr="",
        )

    result = Microsolvator.run(
        solute=solute,
        solvent=solvent,
        config=config,
        constrained_indices=[1, 2],
        constrain_solute=True,
        keep_temps=True,
        working_directory=tmp_path,
        run_command=fake_runner,
    )

    assert (tmp_path / ".xcontrol").exists()
    assert "--nopreopt" in result.command
    assert result.best_structure is not None
    assert len(result.ensemble) == 2
    assert result.population_path is not None
    result.ensure_outputs()
