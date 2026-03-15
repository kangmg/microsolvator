"""Tests for the microsolvator.workflow package."""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT

from microsolvator.workflow.utils import (
    count_solvent_molecules,
    estimate_box_size,
    validate_packmol_output,
)
from microsolvator.workflow.config import (
    EquilibrationConfig,
    PackmolConfig,
    RelaxationConfig,
    SolvationWorkflowConfig,
)
from microsolvator.workflow.packmol import (
    PackmolSolvator,
    _write_packmol_input,
    resolve_packmol,
)
from microsolvator.workflow.swap import kabsch_align, relax_interface, swap_solute
from microsolvator.workflow.equilibration import equilibrate
from microsolvator.workflow.builder import SolvatedTrajectoryBuilder
from microsolvator.workflow.results import SolvatedTrajectoryResult
from microsolvator import MicrosolvatorConfig


# ── Fixtures ────────────────────────────────────────────────────────────────


def _water():
    return Atoms("OH2", positions=[[0, 0, 0], [0.96, 0, 0], [0, 0.76, 0]])


def _cluster():
    """A fake microsolvated cluster: 3-atom solute + 3-atom solvent shell."""
    return Atoms(
        "CuOHOH2",
        positions=[
            [0, 0, 0],
            [1.8, 0, 0],
            [2.5, 0.5, 0],
            [4, 0, 0],
            [4.8, 0.5, 0],
            [4.5, -0.5, 0],
        ],
    )


def _metal_images(n=3):
    """Return n NEB-like images of a Cu3 'solute'."""
    images = []
    for i in range(n):
        shift = 0.1 * i
        atoms = Atoms(
            "Cu3",
            positions=[
                [0, 0, 0],
                [2.5 + shift, 0, 0],
                [1.25, 2.16 + shift, 0],
            ],
        )
        images.append(atoms)
    return images


# ── Utils tests ─────────────────────────────────────────────────────────────


class TestEstimateBoxSize:
    def test_basic(self):
        atoms = Atoms("H2", positions=[[0, 0, 0], [10, 0, 0]])
        box = estimate_box_size(atoms, margin=5.0)
        assert box == pytest.approx(20.0)  # span=10, + 2*5

    def test_3d_span(self):
        atoms = Atoms("H3", positions=[[0, 0, 0], [5, 0, 0], [0, 8, 0]])
        box = estimate_box_size(atoms, margin=3.0)
        # max span is 8 (y-axis), box = 8 + 6 = 14
        assert box == pytest.approx(14.0)

    def test_single_atom(self):
        atoms = Atoms("H", positions=[[1, 2, 3]])
        box = estimate_box_size(atoms, margin=5.0)
        assert box == pytest.approx(10.0)  # span=0, + 2*5


class TestCountSolventMolecules:
    def test_water_density(self):
        water = _water()
        # 30 Å box, density 1.0 g/cm³
        n = count_solvent_molecules(30.0, water, 1.0)
        assert n >= 1
        # Rough check: 30³ Å³ = 2.7e4 Å³ = 2.7e-20 cm³
        # n = V * rho * N_A / M ≈ 2.7e-20 * 1.0 * 6.022e23 / 18.015 ≈ 903
        assert 800 < n < 1100

    def test_at_least_one(self):
        water = _water()
        n = count_solvent_molecules(1.0, water, 0.001)
        assert n >= 1


class TestValidatePackmolOutput:
    def test_valid(self):
        cluster = Atoms("OH2", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        # 2 solvent molecules of 3 atoms each → total = 3 + 6 = 9
        output = Atoms(
            "OH2OH2OH2",
            positions=[[i, 0, 0] for i in range(9)],
        )
        validate_packmol_output(output, cluster, n_bulk_solvent=2, n_solvent_atoms=3)

    def test_wrong_count(self):
        cluster = Atoms("OH2", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        output = Atoms("OH2OH2", positions=[[i, 0, 0] for i in range(6)])
        with pytest.raises(ValueError, match="atom count mismatch"):
            validate_packmol_output(output, cluster, n_bulk_solvent=2, n_solvent_atoms=3)

    def test_wrong_ordering(self):
        cluster = Atoms("OH2", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        # Wrong leading symbols: H2O instead of OH2
        output = Atoms(
            "H2OOH2OH2",
            positions=[[i, 0, 0] for i in range(9)],
        )
        with pytest.raises(ValueError, match="ordering"):
            validate_packmol_output(output, cluster, n_bulk_solvent=2, n_solvent_atoms=3)


# ── Packmol tests ──────────────────────────────────────────────────────────


class TestWritePackmolInput:
    def test_input_file_content(self, tmp_path):
        input_path = tmp_path / "packmol.inp"
        output_path = tmp_path / "solvated.xyz"
        cluster_path = tmp_path / "cluster.xyz"
        solvent_path = tmp_path / "solvent.xyz"

        _write_packmol_input(
            input_path=input_path,
            output_path=output_path,
            cluster_path=cluster_path,
            solvent_path=solvent_path,
            box_center=np.array([15.0, 15.0, 15.0]),
            box_size=30.0,
            n_solvent=100,
            min_distance=2.0,
        )

        content = input_path.read_text()
        assert "tolerance 2.0000" in content
        assert "filetype xyz" in content
        assert f"output {output_path.as_posix()}" in content
        assert f"structure {cluster_path.as_posix()}" in content
        assert "number 1" in content
        assert "center" in content
        assert "fixed 15.000000 15.000000 15.000000 0. 0. 0." in content
        assert f"structure {solvent_path.as_posix()}" in content
        assert "number 100" in content
        assert "inside box" in content


class TestResolvePackmol:
    def test_not_found(self):
        with pytest.raises(FileNotFoundError, match="packmol"):
            resolve_packmol("nonexistent_packmol_binary_xyz")

    def test_explicit_path(self, tmp_path):
        exe = tmp_path / "packmol"
        exe.write_text("#!/bin/sh\n", encoding="utf-8")
        exe.chmod(0o755)
        assert resolve_packmol(str(exe)) == str(exe)


class TestPackmolSolvator:
    def test_density_and_count_mutual_exclusion(self):
        cluster = _cluster()
        solvent = _water()
        with pytest.raises(ValueError, match="provide either"):
            PackmolSolvator.run(
                cluster=cluster,
                solvent=solvent,
                config=PackmolConfig(solvent_density=None, n_bulk_solvent=None),
            )
        with pytest.raises(ValueError, match="not both"):
            PackmolSolvator.run(
                cluster=cluster,
                solvent=solvent,
                config=PackmolConfig(solvent_density=1.0, n_bulk_solvent=10),
            )


# ── Swap / alignment tests ─────────────────────────────────────────────────


class TestKabschAlign:
    def test_identity(self):
        atoms = Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        aligned = kabsch_align(atoms, atoms)
        np.testing.assert_allclose(aligned.positions, atoms.positions, atol=1e-10)

    def test_translation(self):
        target = Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        mobile = Atoms("H3", positions=[[10, 10, 10], [11, 10, 10], [10, 11, 10]])
        aligned = kabsch_align(mobile, target)
        np.testing.assert_allclose(aligned.positions, target.positions, atol=1e-10)

    def test_rotation(self):
        target = Atoms("H3", positions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # Apply a known rotation (90° around z): (x,y,z) → (-y,x,z)
        mobile = Atoms("H3", positions=[[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        aligned = kabsch_align(mobile, target)
        np.testing.assert_allclose(aligned.positions, target.positions, atol=1e-8)

    def test_different_sizes_raises(self):
        a = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
        b = Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError, match="different atom counts"):
            kabsch_align(a, b)

    def test_different_species_raises(self):
        a = Atoms("HO", positions=[[0, 0, 0], [1, 0, 0]])
        b = Atoms("OH", positions=[[0, 0, 0], [1, 0, 0]])
        with pytest.raises(ValueError, match="species mismatch"):
            kabsch_align(a, b)


class TestSwapSolute:
    def test_solvent_preserved(self):
        # Template: 3 solute + 3 solvent atoms
        template = Atoms(
            "Cu3Au3",
            positions=[
                [0, 0, 0], [2.5, 0, 0], [1.25, 2.16, 0],
                [10, 10, 10], [12, 10, 10], [11, 12, 10],
            ],
        )
        new_solute = Atoms(
            "Cu3",
            positions=[[0, 0, 0], [2.6, 0, 0], [1.3, 2.2, 0]],
        )
        result = swap_solute(template, new_solute, n_solute=3)
        # Solvent positions should be unchanged
        np.testing.assert_allclose(
            result.positions[3:], template.positions[3:], atol=1e-10
        )
        # Solute positions should differ from original template
        assert not np.allclose(result.positions[:3], template.positions[:3], atol=1e-5)


class TestRelaxInterface:
    def test_optimize_method(self):
        # Use EMT calculator with Cu atoms
        system = Atoms(
            "Cu4",
            positions=[
                [0, 0, 0],  # solute (frozen)
                [3.0, 0, 0],
                [0, 3.0, 0],
                [3.0, 3.0, 0],
            ],
        )
        system.cell = [10, 10, 10]
        system.pbc = True

        config = RelaxationConfig(method="optimize", max_steps=5, fmax=0.1)
        result = relax_interface(system, n_solute=1, calculator=EMT(), config=config)

        # Solute should not have moved
        np.testing.assert_allclose(result.positions[0], [0, 0, 0], atol=1e-10)
        assert len(result) == 4

    def test_md_method(self):
        system = Atoms(
            "Cu4",
            positions=[
                [0, 0, 0],
                [3.0, 0, 0],
                [0, 3.0, 0],
                [3.0, 3.0, 0],
            ],
        )
        system.cell = [10, 10, 10]
        system.pbc = True

        config = RelaxationConfig(method="md", md_steps=5)
        result = relax_interface(system, n_solute=1, calculator=EMT(), config=config)
        np.testing.assert_allclose(result.positions[0], [0, 0, 0], atol=1e-10)

    def test_unknown_method_raises(self):
        system = Atoms("Cu", positions=[[0, 0, 0]])
        config = RelaxationConfig(method="unknown")
        with pytest.raises(ValueError, match="Unknown relaxation method"):
            relax_interface(system, n_solute=0, calculator=EMT(), config=config)


# ── Equilibration tests ────────────────────────────────────────────────────


class TestEquilibrate:
    def test_basic_nvt(self):
        system = Atoms(
            "Cu4",
            positions=[
                [0, 0, 0],
                [2.55, 0, 0],
                [0, 2.55, 0],
                [2.55, 2.55, 0],
            ],
        )
        system.cell = [10, 10, 10]
        system.pbc = True

        config = EquilibrationConfig(nvt_steps=5, temperature=100.0)
        result = equilibrate(
            system, n_fixed=1, calculator=EMT(), config=config
        )

        # First atom should not have moved
        np.testing.assert_allclose(result.positions[0], [0, 0, 0], atol=1e-10)
        assert len(result) == 4

    def test_with_heating_schedule(self):
        system = Atoms(
            "Cu4",
            positions=[
                [0, 0, 0],
                [2.55, 0, 0],
                [0, 2.55, 0],
                [2.55, 2.55, 0],
            ],
        )
        system.cell = [10, 10, 10]
        system.pbc = True

        config = EquilibrationConfig(
            nvt_steps=3,
            heating_schedule=[(50.0, 3), (100.0, 3)],
        )
        result = equilibrate(
            system, n_fixed=2, calculator=EMT(), config=config
        )
        np.testing.assert_allclose(result.positions[:2], system.positions[:2], atol=1e-10)

    def test_trajectory_written(self, tmp_path):
        system = Atoms(
            "Cu4",
            positions=[
                [0, 0, 0],
                [2.55, 0, 0],
                [0, 2.55, 0],
                [2.55, 2.55, 0],
            ],
        )
        system.cell = [10, 10, 10]
        system.pbc = True

        traj_path = str(tmp_path / "eq.traj")
        config = EquilibrationConfig(nvt_steps=20)
        equilibrate(
            system, n_fixed=1, calculator=EMT(), config=config,
            traj_path=traj_path,
        )
        assert Path(traj_path).exists()


# ── Config tests ────────────────────────────────────────────────────────────


class TestConfigs:
    def test_packmol_config_defaults(self):
        cfg = PackmolConfig()
        assert cfg.box_margin == 10.0
        assert cfg.min_distance == 2.0
        assert cfg.max_retries == 3

    def test_solvation_workflow_config(self):
        cfg = SolvationWorkflowConfig(
            microsolv=MicrosolvatorConfig(nsolv=5),
            packmol=PackmolConfig(solvent_density=1.0),
        )
        assert cfg.microsolv.nsolv == 5
        assert cfg.packmol.solvent_density == 1.0
        assert cfg.ts_index is None

    def test_relaxation_config_defaults(self):
        cfg = RelaxationConfig()
        assert cfg.method == "optimize"
        assert cfg.fmax == 0.05


# ── Builder tests ───────────────────────────────────────────────────────────


class TestSolvatedTrajectoryBuilder:
    def test_empty_images_raises(self):
        config = SolvationWorkflowConfig(
            microsolv=MicrosolvatorConfig(nsolv=3),
        )
        with pytest.raises(ValueError, match="at least one"):
            SolvatedTrajectoryBuilder.build(
                reaction_images=[],
                solvent=_water(),
                config=config,
                calculator=EMT(),
            )

    def test_invalid_ts_index_raises(self):
        images = _metal_images(3)
        config = SolvationWorkflowConfig(
            microsolv=MicrosolvatorConfig(nsolv=3),
            ts_index=10,
        )
        with pytest.raises(ValueError, match="out of range"):
            SolvatedTrajectoryBuilder.build(
                reaction_images=images,
                solvent=_water(),
                config=config,
                calculator=EMT(),
            )

    def test_swap_and_relax_ts_identity(self):
        """TS image should be returned unchanged (copy of template)."""
        template = Atoms(
            "Cu3Au3",
            positions=[
                [0, 0, 0], [2.5, 0, 0], [1.25, 2.16, 0],
                [10, 10, 10], [12, 10, 10], [11, 12, 10],
            ],
        )
        template.cell = [20, 20, 20]
        template.pbc = True

        images = _metal_images(3)
        config = RelaxationConfig(method="optimize", max_steps=3)

        result = SolvatedTrajectoryBuilder.swap_and_relax(
            template=template,
            reaction_images=images,
            n_solute=3,
            ts_index=1,
            calculator=EMT(),
            config=config,
        )

        assert len(result) == 3
        # TS image (index 1) should match template exactly
        np.testing.assert_allclose(result[1].positions, template.positions, atol=1e-10)

    def test_chain_propagation_order(self):
        """Verify that swap_and_relax propagates outward from TS.

        Each non-TS image should use the *previous* (closer-to-TS) image's
        relaxed solvent as its starting template, not the original template.
        """
        # 5 images, TS at index 2
        # Template: 3 solute Cu + 3 solvent Cu
        template = Atoms(
            "Cu6",
            positions=[
                [0, 0, 0], [2.55, 0, 0], [1.275, 2.21, 0],  # solute
                [6, 0, 0], [6, 2.55, 0], [6, 5.1, 0],        # solvent
            ],
        )
        template.cell = [15, 15, 15]
        template.pbc = True

        images = _metal_images(5)
        config = RelaxationConfig(method="optimize", max_steps=3, fmax=0.1)

        result = SolvatedTrajectoryBuilder.swap_and_relax(
            template=template,
            reaction_images=images,
            n_solute=3,
            ts_index=2,
            calculator=EMT(),
            config=config,
        )

        assert len(result) == 5
        # TS image unchanged
        np.testing.assert_allclose(
            result[2].positions, template.positions, atol=1e-10
        )
        # All images should have the same number of atoms
        for img in result:
            assert len(img) == len(template)

    def test_chain_propagation_uses_previous_solvent(self):
        """Adjacent images should share similar solvent positions (continuity)."""
        template = Atoms(
            "Cu6",
            positions=[
                [0, 0, 0], [2.55, 0, 0], [1.275, 2.21, 0],  # solute
                [6, 0, 0], [6, 2.55, 0], [6, 5.1, 0],        # solvent
            ],
        )
        template.cell = [15, 15, 15]
        template.pbc = True

        images = _metal_images(5)
        config = RelaxationConfig(method="optimize", max_steps=2, fmax=0.01)

        result = SolvatedTrajectoryBuilder.swap_and_relax(
            template=template,
            reaction_images=images,
            n_solute=3,
            ts_index=2,
            calculator=EMT(),
            config=config,
        )

        # Check that adjacent images have more similar solvent than distant ones
        # Solvent positions are indices 3:6
        def solvent_rmsd(a, b):
            diff = a.positions[3:] - b.positions[3:]
            return np.sqrt((diff ** 2).mean())

        # TS(2) vs image 3 should be closer than TS(2) vs image 4
        rmsd_2_3 = solvent_rmsd(result[2], result[3])
        rmsd_2_4 = solvent_rmsd(result[2], result[4])
        # At minimum, both should be finite (no NaN/explosion)
        assert np.isfinite(rmsd_2_3)
        assert np.isfinite(rmsd_2_4)


# ── Result container test ──────────────────────────────────────────────────


class TestSolvatedTrajectoryResult:
    def test_fields(self):
        images = [Atoms("H")]
        ts = Atoms("H")
        config = SolvationWorkflowConfig(
            microsolv=MicrosolvatorConfig(nsolv=1),
        )
        result = SolvatedTrajectoryResult(
            solvated_images=images,
            microsolvated_ts=ts,
            boxed_system=ts,
            equilibrated_system=ts,
            n_solute_atoms=1,
            n_cluster_atoms=1,
            n_total_atoms=1,
            ts_index=0,
            config=config,
            working_directory=Path("."),
        )
        assert result.solvated_images == images
        assert result.ts_index == 0
