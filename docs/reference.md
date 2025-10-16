# Microsolvator Full Reference

This reference documents the complete public surface of the `microsolvator` package, covering configuration options, execution semantics, result containers, support helpers, and binary management utilities. For workflow-oriented guidance see `docs/usage.md`.

## Package Layout

| Module | Purpose |
| --- | --- |
| `microsolvator.config` | Configuration objects for CREST microsolvation runs. |
| `microsolvator.runner` | High-level orchestration via `Microsolvator.run`. |
| `microsolvator.command` | Construction of CREST command-line arguments. |
| `microsolvator.results` | Structured run outputs and helpers. |
| `microsolvator.support` | Implicit solvent capability tables and queries. |
| `microsolvator.install` | Binary resolution and installer utilities. |

All primary symbols are re-exported through `microsolvator.__init__` for convenience.

## Configuration (`MicrosolvatorConfig`)

```python
from microsolvator import MicrosolvatorConfig
```

| Field | Type | Default | Description |
| --- | --- | --- | --- |
| `nsolv` | `int` | _(required)_ | Number of solvent molecules to grow via CREST QCG. Must be greater than zero. |
| `method` | `str` | `"gfn2"` | Method string used when validating implicit solvents. |
| `temperature` | `float` | `298.0` | Target temperature passed to CREST (`--T`). |
| `mdtime` | `float` | `50.0` | Meta-dynamics/MD simulation time in ps (`--mdtime`). |
| `charge` | `int` | `0` | Total system charge (`--chrg`). |
| `uhf` | `int` | `0` | Spin multiplicity parameter (`--uhf`). |
| `crest_executable` | `Optional[str]` | `None` | Explicit CREST binary path. Overrides env/package lookup. |
| `xtb_executable` | `Optional[str]` | `None` | Explicit xTB binary path used for `--xnam`. |
| `implicit_model` | `Optional[str]` | `None` | Implicit solvent model flag (`alpb`, `gbsa`, ...). Must be paired with `implicit_solvent`. |
| `implicit_solvent` | `Optional[str]` | `None` | Solvent identifier (e.g. `"h2o"`). Must be paired with `implicit_model`. |
| `ensemble` | `bool` | `True` | Add `--ensemble` when true. |
| `nopreopt` | `bool` | `False` | Adds `--nopreopt`. Automatically enabled when constraints are written. |
| `additional_flags` | `Sequence[str]` | `()` | Extra CREST flags appended verbatim. |

### Helper Methods

- `build_flag_list() -> list[str]`: returns the ordered set of CREST flags derived from the configuration (excluding solute/solvent filenames and resolved binary paths).
- `from_kwargs(**kwargs) -> MicrosolvatorConfig`: convenience alternative to direct instantiation.

## Execution (`Microsolvator.run`)

```python
from microsolvator import Microsolvator

result = Microsolvator.run(
    solute=solute_atoms,
    solvent=solvent_atoms,
    config=config,
    constrained_indices=[...],
    constrain_solute=True,
    keep_temps=False,
    working_directory=None,
    run_command=None,
    prepare_only=False,
)
```

| Parameter | Type | Description |
| --- | --- | --- |
| `solute` | `ase.Atoms` | Solute structure used as CREST starting point. |
| `solvent` | `ase.Atoms` | Solvent molecule to be added repeatedly. |
| `config` | `MicrosolvatorConfig` | Execution configuration (flags, binaries, implicit solvent options, etc.). |
| `constrained_indices` | `Optional[Sequence[int]]` | 1-based atom indices to constrain via `.xcontrol`. Merged with `constrain_solute`. |
| `constrain_solute` | `bool` | When true, every solute atom is constrained. Automatically triggers `--nopreopt` unless already enabled. |
| `keep_temps` | `bool` | Persist temporary directory after run. Automatically enabled when `prepare_only=True` and no `working_directory` is provided. |
| `working_directory` | `Optional[pathlib.Path]` | Directory for input/output files. Created if missing. When provided, caller manages cleanup. |
| `run_command` | `Optional[Callable]` | Optional executor callback. Preferred signature `(command, workdir, env, log_path)` returning `subprocess.CompletedProcess`. Callbacks accepting fewer arguments (e.g. without `env`/`log_path`) remain supported. |
| `prepare_only` | `bool` | When true, inputs are generated and the command is printed/returned without executing CREST. |
| `log_file` | `Optional[str]` | Log file name (relative to working directory) or absolute path. Defaults to `"crest_run.log"`. Set to `None`/empty string to disable streaming logs. |

### Execution Behaviour

1. ASE `Atoms` objects are serialized to `solute.xyz` and `solvent.xyz` (absolute paths) inside the working directory.
2. Constraint indices are merged (from `constrained_indices` and `constrain_solute`). If any exist, `.xcontrol` is written and `--nopreopt` is injected if not already set.
3. CREST and xTB executables are resolved in the following order: explicit config path → environment variable (`CREST_BIN` / `XTB_BIN`) → package `_bin` directory → system `PATH`. Resolved values are converted to absolute paths.
4. `build_crest_command` assembles the final command list using absolute paths for binaries and input files, always appending `--xnam /abs/path/to/xtb`.
5. A subprocess environment is created with `CREST_BIN` and `XTB_BIN` set to the resolved paths.
6. Unless `prepare_only=True`, the command is executed (via either the provided `run_command` or the default executor).
7. During execution, stdout/stderr are streamed to `log_file` (default `crest_run.log`) with `[OUT]`/`[ERR]` prefixes while also being captured in memory for the returned result. Custom `run_command` callbacks receive the resolved log path as the fourth argument and may override logging behaviour.

When `prepare_only=True`, CREST is **not** launched. Inputs remain on disk, the command is printed as a shell string, and the returned result is marked as `executed=False`.

## Results (`MicrosolvationResult`)

Attributes exposed on the returned object:

| Attribute | Type | Description |
| --- | --- | --- |
| `command` | `Sequence[str]` | Tokenized CREST command with absolute paths and `--xnam`. |
| `shell_command` | `str` | Shell-escaped version of `command` (`shlex.join`). |
| `working_directory` | `pathlib.Path` | Directory containing input/output files. |
| `best_structure` | `Optional[ase.Atoms]` | Parsed `crest_best.xyz`, if available. |
| `ensemble` | `list[ase.Atoms]` | Parsed `full_ensemble.xyz` structures. Empty when not produced or `prepare_only=True`. |
| `final` | `Optional[ase.Atoms]` | Parsed `grow/cluster.xyz` (final cluster snapshot), when available. |
| `traj` | `list[ase.Atoms]` | Parsed `grow/qcg_grow.xyz` trajectory frames (may be empty). |
| `population_path` | `Optional[pathlib.Path]` | Path to `full_population.dat` when present. |
| `stdout` / `stderr` | `str` | Captured CREST output streams (empty in prepare-only mode). |
| `executed` | `bool` | Indicates whether CREST was actually run. |

Methods:

- `ensure_outputs()`: Raises `ValueError` if `executed=True` and no structures were parsed. No-op in prepare-only mode.

## Implicit Solvent Helpers (`microsolvator.support`)

- `supports_implicit_solvent(*, method: str, model: str, solvent: str) -> bool`
  - Returns whether the combination is supported according to packaged tables.

- `list_supported_implicit_solvents(*, method: Optional[str] = None, model: Optional[str] = None) -> dict`
  - Provides a nested dictionary of supported solvents filtered by method/model (if given).

- `validate_implicit_choice(*, method: str, model: str, solvent: str) -> None`
  - Raises `ValueError` when the combination is unsupported. Invoked automatically by `Microsolvator.run` when implicit solvation is requested.

## Binary Utilities (`microsolvator.install`)

| Symbol | Description |
| --- | --- |
| `install_crest(url=DEFAULT_CREST_URL, force=False) -> Path` | Downloads the official CREST tarball, extracts the binary into `microsolvator/_bin`, removes archives, and ensures executable bits. |
| `install_xtb(url=DEFAULT_XTB_URL, force=False) -> Path` | Same as above for xTB. |
| `resolve_crest_binary(explicit_path: Optional[str]) -> str` | Resolves CREST executable (explicit path → env → package `_bin` → `PATH`). Returns absolute path when found. |
| `resolve_xtb_binary(explicit_path: Optional[str]) -> str` | xTB analogue of the resolver above. |
| `PACKAGE_BIN_DIR: Path` | Directory where bundled binaries are stored (created on demand). |

Both installers delete temporary `.tar.xz` archives after extraction. Use `force=True` to overwrite existing binaries.

## Command Generation Example

```python
from pathlib import Path
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

solute = read("benzoic_acid.xyz")
solvent = read("water.xyz")

config = MicrosolvatorConfig(
    nsolv=3,
    implicit_model="alpb",
    implicit_solvent="h2o",
    crest_executable="/content/crest/crest",
    xtb_executable="/content/xtb-dist/bin/xtb",
)

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    working_directory=Path("/content"),
    constrain_solute=True,
    prepare_only=True,
)

print(result.shell_command)
```

Possible output:

```
/content/crest/crest /content/solute.xyz --qcg /content/solvent.xyz --nsolv 3 --T 298.0 --mdtime 50.0 --ensemble --alpb h2o --chrg 0 --uhf 0 --nopreopt --xnam /content/xtb-dist/bin/xtb
```

The command string can be copied directly into a shell or scheduler script.

## Prepare-Only Workflow Summary

1. Set `prepare_only=True` (optionally alongside `working_directory`).
2. Inspect `result.shell_command` and generated input files (`solute.xyz`, `solvent.xyz`, `.xcontrol`).
3. Execute the printed command manually when ready.

This mode is useful for dry runs, integration with custom schedulers, or debugging CREST inputs without launching the actual computation.
