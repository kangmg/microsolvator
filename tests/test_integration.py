"""Integration tests that require real CREST and xTB binaries.

These tests are skipped unless the MICROSOLVATOR_INTEGRATION environment
variable is set (e.g. in CI).  They exercise install_crest / install_xtb,
binary resolution, and a real CREST QCG run with a tiny molecule.
"""

from __future__ import annotations

import os
import subprocess

import pytest
from ase import Atoms

from microsolvator import (
    Microsolvator,
    MicrosolvatorConfig,
    install_crest,
    install_xtb,
    resolve_crest_binary,
    resolve_xtb_binary,
)

INTEGRATION = os.environ.get("MICROSOLVATOR_INTEGRATION", "")
pytestmark = pytest.mark.skipif(
    not INTEGRATION,
    reason="Set MICROSOLVATOR_INTEGRATION=1 to run integration tests",
)


# ------------------------------------------------------------------
# 1. install_crest / install_xtb
# ------------------------------------------------------------------

class TestInstallBinaries:

    def test_install_xtb(self):
        path = install_xtb(force=True)
        assert path.exists()
        assert path.name == "xtb"
        result = subprocess.run(
            [str(path), "--version"],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0 or "xtb" in (result.stdout + result.stderr).lower()

    def test_install_crest(self):
        path = install_crest(force=True)
        assert path.exists()
        assert path.name == "crest"
        result = subprocess.run(
            [str(path), "--version"],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0 or "crest" in (result.stdout + result.stderr).lower()

    def test_reinstall_without_force_returns_existing(self):
        first = install_xtb(force=True)
        second = install_xtb(force=False)
        assert first == second


# ------------------------------------------------------------------
# 2. resolve_*_binary finds installed binaries
# ------------------------------------------------------------------

class TestResolveBinaries:

    def test_resolve_crest_finds_installed(self):
        install_crest(force=True)
        path = resolve_crest_binary(None)
        assert "crest" in path

    def test_resolve_xtb_finds_installed(self):
        install_xtb(force=True)
        path = resolve_xtb_binary(None)
        assert "xtb" in path


# ------------------------------------------------------------------
# 3. Real CREST QCG run (tiny molecule, fast)
# ------------------------------------------------------------------

class TestCrestQcgRun:

    def test_h2o_microsolvation_1_solvent(self, tmp_path):
        """Run a real QCG grow with 1 solvent molecule — fast sanity check."""
        install_crest(force=False)
        install_xtb(force=False)

        solute = Atoms(
            "OH2",
            positions=[[0.0, 0.0, 0.1173], [-0.7572, 0.0, -0.4692], [0.7572, 0.0, -0.4692]],
        )
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method="gfn2",
            threads=1,
            mdtime=5.0,        # very short for CI speed
            ensemble=False,    # skip ensemble search — just grow
        )

        result = Microsolvator.run(
            solute=solute,
            solvent=solvent,
            config=config,
            working_directory=tmp_path,
            keep_temps=True,
        )

        # grow must produce a final cluster
        assert result.final is not None, "grow/cluster.xyz was not produced"
        assert len(result.final) > len(solute), "final cluster should be larger than solute"

        # stdout should contain CREST output
        assert result.stdout, "CREST produced no stdout"

    def test_h2o_microsolvation_with_ensemble(self, tmp_path):
        """Run grow + ensemble with 1 solvent molecule."""
        install_crest(force=False)
        install_xtb(force=False)

        solute = Atoms(
            "OH2",
            positions=[[0.0, 0.0, 0.1173], [-0.7572, 0.0, -0.4692], [0.7572, 0.0, -0.4692]],
        )
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method="gfn2",
            threads=1,
            mdtime=5.0,
            ensemble=True,
        )

        result = Microsolvator.run(
            solute=solute,
            solvent=solvent,
            config=config,
            working_directory=tmp_path,
            keep_temps=True,
        )

        assert result.final is not None
        assert result.best_structure is not None, "crest_best.xyz was not produced"
        assert len(result.ensemble) >= 1, "ensemble should have at least 1 conformer"

    def test_h2o_with_constraints(self, tmp_path):
        """Run grow with solute constraints."""
        install_crest(force=False)
        install_xtb(force=False)

        solute = Atoms(
            "OH2",
            positions=[[0.0, 0.0, 0.1173], [-0.7572, 0.0, -0.4692], [0.7572, 0.0, -0.4692]],
        )
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method="gfn2",
            threads=1,
            mdtime=5.0,
            ensemble=False,
        )

        result = Microsolvator.run(
            solute=solute,
            solvent=solvent,
            config=config,
            constrain_solute=True,
            working_directory=tmp_path,
            keep_temps=True,
        )

        assert result.final is not None
        assert "--nopreopt" in result.command
        assert (tmp_path / ".xcontrol").exists()

    def test_prepare_only_does_not_run(self, tmp_path):
        """prepare_only should write inputs but not execute CREST."""
        install_crest(force=False)
        install_xtb(force=False)

        solute = Atoms(
            "OH2",
            positions=[[0.0, 0.0, 0.1173], [-0.7572, 0.0, -0.4692], [0.7572, 0.0, -0.4692]],
        )
        solvent = solute.copy()

        config = MicrosolvatorConfig(nsolv=1, method="gfn2")

        result = Microsolvator.run(
            solute=solute,
            solvent=solvent,
            config=config,
            working_directory=tmp_path,
            prepare_only=True,
        )

        assert result.executed is False
        assert (tmp_path / "solute.xyz").exists()
        assert (tmp_path / "solvent.xyz").exists()
        assert result.best_structure is None
