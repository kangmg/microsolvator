# Microsolvator Usage Guide

This guide clarifies how `Microsolvator` converts user input into the CREST command line and how to customize execution.

## Command Construction

`Microsolvator.run` assembles a list of shell tokens that will be passed to `subprocess.run`. The final list is accessible via `MicrosolvationResult.command`.

Example:

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

solute = read("benzoic_acid.xyz")
solvent = read("water.xyz")

config = MicrosolvatorConfig(
    nsolv=3,
    temperature=298.0,
    mdtime=50.0,
    implicit_model="alpb",
    implicit_solvent="h2o",
)

result = Microsolvator.run(solute=solute, solvent=solvent, config=config)
print(result.command)
```

Possible output:

```
['/abs/path/crest', '/abs/path/solute.xyz', '--qcg', '/abs/path/solvent.xyz', '--nsolv', '3', '--temp', '298', '--T', '1', '--mdtime', '50.0', '--enslvl', 'gfn2', '--ensemble', '--alpb', 'h2o', '--chrg', '0', '--uhf', '0', '--xnam', '/abs/path/xtb']
```

The command always begins with an auto-resolved CREST executable. The resolution order is:

1. `config.crest_executable` if explicitly provided.
2. Environment variable `CREST_BIN` if it points to an existing file.
3. A binary stored under `microsolvator/_bin/` (populated via `install_crest`).
4. Fallback to `"crest"`, allowing the system `PATH` to resolve it.

xTB paths follow the same order using `config.xtb_executable`, `XTB_BIN`, the package directory (via `install_xtb`), and finally `"xtb"`. The resolved path is injected both via `--xnam` (so CREST sees it explicitly) and the `XTB_BIN` environment variable.

Threads (CREST `--T`) default to `config.threads` (1 if unspecified). Temperature is mapped to `--temp` instead of `--T`.

## Custom Execution (`run_command`)

`Microsolvator.run` accepts a `run_command` callback if you need to redirect command execution (e.g., mocking in tests). The preferred signature is:

```python
def run_command(
    command: Sequence[str],
    workdir: Path,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    ...
```

- `command` is the list described above.
- `workdir` is the directory containing the generated input files (solute.xyz, solvent.xyz, optional `.xcontrol`).
- `env` contains the environment passed to `subprocess.run`, pre-populated with `CREST_BIN` and `XTB_BIN` based on the resolution logic.
- The callback must return a `subprocess.CompletedProcess` with `stdout`, `stderr`, and `returncode` set.

When `run_command` is omitted, `subprocess.run` is used with `check=True`, `capture_output=True`, `text=True`, and the environment described above.

Callbacks that only accept `(command, workdir)` also work; the environment argument is optional.

### Live Logging

By default, CREST stdout/stderr are streamed to `crest_run.log` inside the working directory while the process runs. Supply `log_file="custom.log"` (or an absolute path) to change the destination, or `log_file=None`/`""` to disable streaming.

Custom `run_command` callbacks receive the resolved log path as the fourth argument and may implement bespoke logging.

## Constraint Handling

If `constrain_solute=True` or `constrained_indices` are provided, an `.xcontrol` file is written in the working directory and `--nopreopt` is appended to the command (unless already enabled in the config).

## Working Directory Control

- `working_directory`: points microsolvator to an existing folder. Temporary files persist in this location.
- `keep_temps=True`: generates a unique temporary directory and keeps it after the run completes.
- Default (no `working_directory`, `keep_temps=False`): temporary files are discarded automatically.

The working directory always contains the files passed to CREST (`solute.xyz`, `solvent.xyz`, `.xcontrol` if any). Output files (`crest_best.xyz`, `full_ensemble.xyz`, `full_population.dat`) are read from the same directory.

## Managing Bundled Binaries

The helper functions `microsolvator.install.install_crest` and `microsolvator.install.install_xtb` download the official Linux tarballs, extract the executables into `microsolvator/_bin/`, ensure they are executable, and remove the downloaded `.tar.xz` archives. Once installed, the resolvers automatically pick them up without additional configuration.

## Result Introspection

`MicrosolvationResult` provides:

- `command`: list of shell tokens used to invoke CREST.
- `shell_command`: shell-escaped string representation of the command.
- `best_structure`: lowest-energy structure as an ASE `Atoms` instance, if produced.
- `ensemble`: list of `Atoms` for the entire ensemble (from `full_ensemble.xyz`).
- `final`: final cluster structure from `grow/cluster.xyz`, when available.
- `traj`: list of trajectory frames parsed from `grow/qcg_grow.xyz`.
- `population_path`: Path to `full_population.dat` when present.
- `stdout` / `stderr`: captured CREST output strings.
- `executed`: boolean flag indicating whether CREST was run (False in `prepare_only` mode).

Use `result.ensure_outputs()` to assert that at least one structure was parsed.

## Prepare-Only Mode

Set `prepare_only=True` to generate input files (`solute.xyz`, `solvent.xyz`, optional `.xcontrol`) and print the fully qualified CREST command without executing it:

```python
result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrain_solute=True,
    working_directory=Path("./scratch"),
    prepare_only=True,
)

print(result.shell_command)
```

When no explicit working directory is given, `prepare_only=True` implicitly enables `keep_temps=True` so that the generated files remain on disk for inspection. `MicrosolvatorResult.executed` will be `False`, and `stdout` / `stderr` are empty strings.
