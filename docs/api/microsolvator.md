# microsolvator

## `MicrosolvatorConfig`

```python
from microsolvator import MicrosolvatorConfig
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `nsolv` | `int` | *(required)* | Number of solvent molecules to add (CREST `--nsolv`) |
| `method` | `str` | `"gfn2"` | Level of theory: `"gfn0"`, `"gfn1"`, `"gfn2"`, `"gfnff"` |
| `temperature` | `float` | `298.0` | Temperature in K (`--temp`) |
| `mdtime` | `float` | `50.0` | MD simulation time in ps (`--mdtime`) |
| `threads` | `int` | `1` | Number of threads (`--T`) |
| `charge` | `int` | `0` | Total charge (`--chrg`) |
| `uhf` | `int` | `0` | Unpaired electrons (`--uhf`) |
| `crest_executable` | `str \| None` | `None` | Explicit CREST binary path |
| `xtb_executable` | `str \| None` | `None` | Explicit xTB binary path |
| `implicit_model` | `str \| None` | `None` | Implicit solvent model: `"alpb"`, `"gbsa"` |
| `implicit_solvent` | `str \| None` | `None` | Solvent name, e.g. `"h2o"` |
| `ensemble` | `bool` | `True` | Run ensemble search after grow (`--ensemble`) |
| `nopreopt` | `bool` | `False` | Skip pre-optimisation (`--nopreopt`) |
| `additional_flags` | `Sequence[str]` | `()` | Extra flags appended verbatim to the CREST command |

---

## `Microsolvator.run`

```python
Microsolvator.run(
    *,
    solute: Atoms,
    solvent: Atoms,
    config: MicrosolvatorConfig,
    constrained_indices: Sequence[int] | None = None,
    constrain_solute: bool = False,
    keep_temps: bool = False,
    working_directory: Path | None = None,
    run_command: Callable | None = None,
    prepare_only: bool = False,
    log_file: str | None = "crest_run.log",
) -> MicrosolvationResult
```

| Parameter | Description |
|-----------|-------------|
| `solute` | Solute `Atoms` object |
| `solvent` | Solvent molecule `Atoms` object |
| `config` | `MicrosolvatorConfig` |
| `constrained_indices` | 1-based atom indices to freeze (`.xcontrol`) |
| `constrain_solute` | Freeze all solute atoms |
| `keep_temps` | Keep the auto-generated temp directory |
| `working_directory` | Use this directory instead of a temp dir |
| `run_command` | Custom executor callback |
| `prepare_only` | Generate inputs and print command; do not run CREST |
| `log_file` | Log path (relative to workdir). `None` disables logging |

---

## `MicrosolvationResult`

```python
from microsolvator import MicrosolvationResult
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `best_structure` | `Atoms \| None` | `crest_best.xyz` — lowest-energy conformer (`ensemble=True`) |
| `final` | `Atoms \| None` | `grow/cluster.xyz` — final grow cluster |
| `ensemble` | `list[Atoms]` | `full_ensemble.xyz` — all conformers |
| `traj` | `list[Atoms]` | `grow/qcg_grow.xyz` — grow trajectory frames |
| `population_path` | `Path \| None` | Path to `full_population.dat` |
| `command` | `Sequence[str]` | Tokenised CREST command |
| `shell_command` | `str` | Shell-escaped command string |
| `working_directory` | `Path` | Directory containing all files |
| `stdout` / `stderr` | `str` | Captured output |
| `executed` | `bool` | `False` in `prepare_only` mode |

**Method:** `ensure_outputs()` — raises `ValueError` if `executed=True` and no structures were parsed.

---

## Implicit solvent helpers

```python
from microsolvator import supports_implicit_solvent, list_supported_implicit_solvents

# Check a specific combination
supports_implicit_solvent(method="gfn2", model="alpb", solvent="h2o")  # bool

# List all supported solvents
list_supported_implicit_solvents(method="gfn2", model="alpb")
```

---

## Binary utilities

```python
from microsolvator import install_crest, install_xtb, resolve_crest_binary, resolve_xtb_binary

install_crest()          # download & install to microsolvator/_bin/
install_xtb()

resolve_crest_binary(None)   # returns resolved absolute path
resolve_xtb_binary("/my/xtb")
```

Resolution order: explicit path → `CREST_BIN`/`XTB_BIN` env var → `microsolvator/_bin/` → system PATH.
