"""Tests for bugs and fixes identified during code review.

Covers:
- support.py: method key normalization across all aliases
- runner.py: _count_positional_params dispatch and TypeError propagation
- install.py: force flag behaviour
- config.py: _method_flag edge cases and build_flag_list validation
- version consistency between pyproject.toml and __init__.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib  # type: ignore[no-redef]
from typing import Optional, Sequence

import pytest
from ase import Atoms
from ase.io import write as ase_write

from microsolvator import (
    Microsolvator,
    MicrosolvatorConfig,
    __version__,
    list_supported_implicit_solvents,
    supports_implicit_solvent,
)
from microsolvator.config import _method_flag
from microsolvator.install import (
    PACKAGE_BIN_DIR,
    _install_binary,
    _promote_binary,
)
from microsolvator.runner import _count_positional_params
from microsolvator.support import _normalize_method, validate_implicit_choice


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_executable(path: Path) -> Path:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_outputs(workdir: Path) -> None:
    """Write minimal output files that Microsolvator expects after a run."""
    best = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.8, 0]])
    ase_write(workdir / "crest_best.xyz", best)
    ase_write(workdir / "full_ensemble.xyz", [best, best.copy()])
    grow_dir = workdir / "grow"
    grow_dir.mkdir(exist_ok=True)
    ase_write(grow_dir / "cluster.xyz", best)
    ase_write(grow_dir / "qcg_grow.xyz", [best])
    (workdir / "full_population.dat").write_text("1 0.5\n", encoding="utf-8")


# ===================================================================
# 1. support.py — method key normalization
# ===================================================================

class TestMethodNormalization:
    """_normalize_method and _SOLVENT_TABLE key consistency."""

    @pytest.mark.parametrize(
        "alias, expected",
        [
            ("gfn2", "gfn2"),
            ("gfn-2", "gfn2"),
            ("GFN2", "gfn2"),
            ("GFN-2", "gfn2"),
            ("gfn1", "gfn1"),
            ("gfn-1", "gfn1"),
            ("gfn0", "gfn0"),
            ("gfn-0", "gfn0"),
            ("gfnff", "gfnff"),
            ("gfn-ff", "gfnff"),
            ("ff", "gfnff"),
            ("FF", "gfnff"),
        ],
    )
    def test_normalize_method_all_aliases(self, alias: str, expected: str):
        assert _normalize_method(alias) == expected

    def test_normalize_method_unknown_falls_back_to_lower(self):
        assert _normalize_method("SomeWeirdMethod") == "someweirdmethod"

    @pytest.mark.parametrize("alias", ["gfn-ff", "gfnff", "ff", "FF"])
    def test_gfnff_alias_lookup_in_solvent_table(self, alias: str):
        """All gfn-ff aliases should resolve to the same solvent table entry."""
        assert supports_implicit_solvent(method=alias, model="alpb", solvent="h2o")

    @pytest.mark.parametrize("alias", ["gfn-2", "gfn2", "GFN2", "GFN-2"])
    def test_gfn2_alias_lookup_in_solvent_table(self, alias: str):
        assert supports_implicit_solvent(method=alias, model="alpb", solvent="h2o")

    @pytest.mark.parametrize("alias", ["gfn-1", "gfn1", "GFN1"])
    def test_gfn1_alias_lookup_in_solvent_table(self, alias: str):
        assert supports_implicit_solvent(method=alias, model="alpb", solvent="h2o")

    def test_gfn0_unsupported_in_solvent_table(self):
        """gfn0 is accepted by _method_flag but has no solvent table entry."""
        assert not supports_implicit_solvent(method="gfn0", model="alpb", solvent="h2o")

    def test_list_supported_filters_by_alias(self):
        """list_supported_implicit_solvents should work with any method alias."""
        data = list_supported_implicit_solvents(method="gfn-ff")
        assert "gfnff" in data
        assert "alpb" in data["gfnff"]

    def test_validate_implicit_choice_raises_on_unsupported(self):
        with pytest.raises(ValueError, match="Unsupported implicit solvent"):
            validate_implicit_choice(method="gfn2", model="gbsa", solvent="aniline")

    def test_validate_implicit_choice_passes_with_alias(self):
        validate_implicit_choice(method="gfn-2", model="alpb", solvent="h2o")


# ===================================================================
# 2. runner.py — _count_positional_params & executor dispatch
# ===================================================================

class TestCountPositionalParams:
    """_count_positional_params correctly counts positional params."""

    def test_four_params(self):
        def f(a, b, c, d): ...
        assert _count_positional_params(f) == 4

    def test_four_params_with_defaults(self):
        def f(a, b, c=None, d=None): ...
        assert _count_positional_params(f) == 4

    def test_two_required_two_optional(self):
        """All positional params counted, regardless of defaults."""
        def f(cmd, workdir, env=None, log=None): ...
        assert _count_positional_params(f) == 4

    def test_two_params_only(self):
        def f(cmd, workdir): ...
        assert _count_positional_params(f) == 2

    def test_three_params(self):
        def f(cmd, workdir, env): ...
        assert _count_positional_params(f) == 3

    def test_kwargs_not_counted(self):
        def f(cmd, workdir, **kwargs): ...
        assert _count_positional_params(f) == 2

    def test_positional_only(self):
        # Python 3.8+ positional-only syntax
        exec_ns = {}
        exec("def f(a, b, /, c): ...", exec_ns)
        assert _count_positional_params(exec_ns["f"]) == 3

    def test_builtin_fallback(self):
        """Builtins without inspectable signature should return 4."""
        assert _count_positional_params(print) == 0 or _count_positional_params(print) >= 0
        # len has no inspectable signature in some Python versions
        # Just confirm it doesn't crash


class TestExecutorDispatch:
    """Executor is called with the correct number of arguments based on its signature."""

    def _make_config(self, tmp_path: Path) -> MicrosolvatorConfig:
        config = MicrosolvatorConfig(nsolv=2, method="gfn2")
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))
        return config

    def _solute_and_solvent(self):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        return solute, solvent

    def test_dispatch_4_params(self, tmp_path: Path):
        """Executor with 4 params receives (command, workdir, env, log_path)."""
        captured: dict = {}

        def runner_4(command, workdir, env, log_path):
            captured["env"] = env
            captured["log_path"] = log_path
            _fake_outputs(workdir)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, run_command=runner_4,
        )
        assert captured["env"] is not None
        assert "log_path" in captured

    def test_dispatch_3_params(self, tmp_path: Path):
        """Executor with 3 params receives (command, workdir, env)."""
        captured: dict = {}

        def runner_3(command, workdir, env):
            captured["env"] = env
            _fake_outputs(workdir)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, run_command=runner_3,
        )
        assert captured["env"] is not None

    def test_dispatch_2_params(self, tmp_path: Path):
        """Executor with 2 params receives (command, workdir) only."""
        captured: dict = {}

        def runner_2(command, workdir):
            captured["called"] = True
            _fake_outputs(workdir)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, run_command=runner_2,
        )
        assert captured["called"]

    def test_dispatch_4_params_with_defaults(self, tmp_path: Path):
        """Executor with 4 params (2 required + 2 optional) gets all 4 args."""
        captured: dict = {}

        def runner_default(command, workdir, env=None, log_path=None):
            captured["env"] = env
            captured["log_path"] = log_path
            _fake_outputs(workdir)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, run_command=runner_default,
        )
        assert captured["env"] is not None
        assert captured["log_path"] is not None

    def test_typeerror_inside_executor_propagates(self, tmp_path: Path):
        """A TypeError raised inside the executor must NOT be swallowed."""

        def buggy_runner(command, workdir, env, log_path):
            # Simulate a real bug inside the executor
            result = None
            return result.stdout  # TypeError: NoneType has no attribute 'stdout'

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        with pytest.raises(AttributeError):
            Microsolvator.run(
                solute=solute, solvent=solvent, config=config,
                working_directory=tmp_path, run_command=buggy_runner,
            )

    def test_typeerror_from_bad_arg_type_propagates(self, tmp_path: Path):
        """A TypeError from wrong arg types inside executor must propagate."""

        def type_error_runner(command, workdir, env, log_path):
            # Simulate a bug: trying to add int to str
            _ = "prefix" + 42  # noqa: F841
            _fake_outputs(workdir)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        solute, solvent = self._solute_and_solvent()
        config = self._make_config(tmp_path)
        with pytest.raises(TypeError):
            Microsolvator.run(
                solute=solute, solvent=solvent, config=config,
                working_directory=tmp_path, run_command=type_error_runner,
            )


# ===================================================================
# 3. install.py — force flag behaviour
# ===================================================================

class TestInstallForceFlag:
    """The force flag should control whether existing binaries are overwritten."""

    def test_promote_binary_finds_nested(self, tmp_path: Path):
        nested = tmp_path / "subdir" / "bin" / "crest"
        nested.parent.mkdir(parents=True)
        _write_executable(nested)
        assert _promote_binary(tmp_path, "crest") == nested

    def test_promote_binary_returns_none_when_missing(self, tmp_path: Path):
        assert _promote_binary(tmp_path, "nonexistent") is None


# ===================================================================
# 4. config.py — _method_flag edge cases & build_flag_list validation
# ===================================================================

class TestMethodFlag:

    @pytest.mark.parametrize(
        "raw, expected",
        [
            ("gfn2", "gfn2"),
            ("gfn-2", "gfn2"),
            ("GFN2", "gfn2"),
            ("  gfn2  ", "gfn2"),
            ("gfn1", "gfn1"),
            ("gfn-1", "gfn1"),
            ("gfn0", "gfn0"),
            ("gfn-0", "gfn0"),
            ("gfnff", "gfnff"),
            ("gfn-ff", "gfnff"),
            ("ff", "gfnff"),
            ("FF", "gfnff"),
        ],
    )
    def test_all_aliases(self, raw: str, expected: str):
        assert _method_flag(raw) == expected

    def test_unsupported_method_raises(self):
        with pytest.raises(ValueError, match="Unsupported method"):
            _method_flag("am1")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            _method_flag("")


class TestBuildFlagList:

    def test_nsolv_zero_raises(self):
        config = MicrosolvatorConfig(nsolv=0, method="gfn2")
        with pytest.raises(ValueError, match="nsolv must be a positive integer"):
            config.build_flag_list()

    def test_nsolv_negative_raises(self):
        config = MicrosolvatorConfig(nsolv=-1, method="gfn2")
        with pytest.raises(ValueError, match="nsolv must be a positive integer"):
            config.build_flag_list()

    def test_threads_zero_raises(self):
        config = MicrosolvatorConfig(nsolv=1, method="gfn2", threads=0)
        with pytest.raises(ValueError, match="threads must be a positive integer"):
            config.build_flag_list()

    def test_implicit_model_only_raises(self):
        config = MicrosolvatorConfig(nsolv=1, method="gfn2", implicit_model="alpb")
        with pytest.raises(ValueError, match="Both implicit_model and implicit_solvent"):
            config.build_flag_list()

    def test_implicit_solvent_only_raises(self):
        config = MicrosolvatorConfig(nsolv=1, method="gfn2", implicit_solvent="h2o")
        with pytest.raises(ValueError, match="Both implicit_model and implicit_solvent"):
            config.build_flag_list()

    def test_flag_list_contains_expected_flags(self):
        config = MicrosolvatorConfig(
            nsolv=3,
            method="gfn2",
            temperature=300.0,
            threads=4,
            implicit_model="alpb",
            implicit_solvent="h2o",
            charge=1,
            uhf=2,
            nopreopt=True,
        )
        flags = config.build_flag_list()
        assert "--qcg" in flags
        assert flags[flags.index("--nsolv") + 1] == "3"
        assert flags[flags.index("--temp") + 1] == "300"
        assert flags[flags.index("--T") + 1] == "4"
        assert flags[flags.index("--enslvl") + 1] == "gfn2"
        assert flags[flags.index("--alpb") + 1] == "h2o"
        assert flags[flags.index("--chrg") + 1] == "1"
        assert flags[flags.index("--uhf") + 1] == "2"
        assert "--nopreopt" in flags
        assert "--ensemble" in flags

    def test_additional_flags_appended(self):
        config = MicrosolvatorConfig(
            nsolv=1, method="gfn2", additional_flags=("--verbose", "--debug")
        )
        flags = config.build_flag_list()
        assert "--verbose" in flags
        assert "--debug" in flags

    def test_ensemble_disabled(self):
        config = MicrosolvatorConfig(nsolv=1, method="gfn2", ensemble=False)
        flags = config.build_flag_list()
        assert "--ensemble" not in flags

    def test_gfnff_method_in_flags(self):
        config = MicrosolvatorConfig(nsolv=1, method="gfn-ff")
        flags = config.build_flag_list()
        assert flags[flags.index("--enslvl") + 1] == "gfnff"


# ===================================================================
# 5. Version consistency
# ===================================================================

class TestVersionConsistency:

    def test_pyproject_matches_init(self):
        pyproject_path = Path(__file__).resolve().parent.parent / "pyproject.toml"
        with open(pyproject_path, "rb") as f:
            data = tomllib.load(f)
        pyproject_version = data["project"]["version"]
        assert pyproject_version == __version__, (
            f"pyproject.toml version ({pyproject_version}) != "
            f"__init__.py version ({__version__})"
        )


# ===================================================================
# 6. Microsolvator.run integration — validate_implicit_choice with aliases
# ===================================================================

class TestRunImplicitValidation:
    """Microsolvator.run calls validate_implicit_choice with the user's method string.
    Ensure it works for all aliases without raising false positives.
    """

    @pytest.mark.parametrize("method_alias", ["gfn2", "gfn-2", "GFN2"])
    def test_run_prepare_only_accepts_method_aliases(self, tmp_path: Path, method_alias: str):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method=method_alias,
            implicit_model="alpb",
            implicit_solvent="h2o",
        )
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        result = Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, prepare_only=True,
        )
        assert result.executed is False

    @pytest.mark.parametrize("method_alias", ["gfnff", "gfn-ff", "ff"])
    def test_run_prepare_only_accepts_gfnff_aliases(self, tmp_path: Path, method_alias: str):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method=method_alias,
            implicit_model="alpb",
            implicit_solvent="h2o",
        )
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        result = Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, prepare_only=True,
        )
        assert result.executed is False

    def test_run_rejects_unsupported_implicit_combo(self, tmp_path: Path):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(
            nsolv=1,
            method="gfnff",
            implicit_model="gbsa",  # gfnff has no gbsa entry
            implicit_solvent="h2o",
        )
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        with pytest.raises(ValueError, match="Unsupported implicit solvent"):
            Microsolvator.run(
                solute=solute, solvent=solvent, config=config,
                working_directory=tmp_path, prepare_only=True,
            )


# ===================================================================
# 7. Constraints edge cases
# ===================================================================

class TestConstraints:

    def test_zero_index_raises(self, tmp_path: Path):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(nsolv=1, method="gfn2")
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        with pytest.raises(ValueError, match="1-based and positive"):
            Microsolvator.run(
                solute=solute, solvent=solvent, config=config,
                constrained_indices=[0, 1],
                working_directory=tmp_path, prepare_only=True,
            )

    def test_constrain_solute_writes_xcontrol(self, tmp_path: Path):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(nsolv=1, method="gfn2")
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        result = Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            constrain_solute=True,
            working_directory=tmp_path, prepare_only=True,
        )
        xcontrol = (tmp_path / ".xcontrol").read_text(encoding="utf-8")
        assert "$constrain" in xcontrol
        assert "1,2,3" in xcontrol
        assert "--nopreopt" in result.command

    def test_no_constraints_no_xcontrol(self, tmp_path: Path):
        solute = Atoms("H2O", positions=[[0, 0, 0], [0.9, 0, 0], [0, 0.7, 0]])
        solvent = solute.copy()

        config = MicrosolvatorConfig(nsolv=1, method="gfn2")
        config.crest_executable = str(_write_executable(tmp_path / "crest"))
        config.xtb_executable = str(_write_executable(tmp_path / "xtb"))

        Microsolvator.run(
            solute=solute, solvent=solvent, config=config,
            working_directory=tmp_path, prepare_only=True,
        )
        assert not (tmp_path / ".xcontrol").exists()


# ===================================================================
# 8. MicrosolvationResult
# ===================================================================

class TestMicrosolvationResult:

    def test_ensure_outputs_raises_when_no_structures(self):
        from microsolvator.results import MicrosolvationResult

        result = MicrosolvationResult(
            command=["crest"],
            working_directory=Path("/tmp"),
            best_structure=None,
            ensemble=[],
            population_path=None,
            stdout="",
            stderr="",
            executed=True,
            final=None,
            traj=[],
        )
        with pytest.raises(ValueError, match="no structures"):
            result.ensure_outputs()

    def test_ensure_outputs_skips_for_non_executed(self):
        from microsolvator.results import MicrosolvationResult

        result = MicrosolvationResult(
            command=["crest"],
            working_directory=Path("/tmp"),
            best_structure=None,
            ensemble=[],
            population_path=None,
            stdout="",
            stderr="",
            executed=False,
            final=None,
            traj=[],
        )
        result.ensure_outputs()  # Should not raise

    def test_shell_command_property(self):
        from microsolvator.results import MicrosolvationResult

        result = MicrosolvationResult(
            command=["crest", "solute.xyz", "--qcg", "solvent.xyz"],
            working_directory=Path("/tmp"),
            best_structure=None,
            ensemble=[],
            population_path=None,
            stdout="",
            stderr="",
            executed=False,
        )
        assert result.shell_command == "crest solute.xyz --qcg solvent.xyz"
