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
['crest', 'solute.xyz', '--qcg', 'solvent.xyz', '--nsolv', '3', '--T', '298.0', '--mdtime', '50.0', '--ensemble', '--alpb', 'h2o', '--chrg', '0', '--uhf', '0']
```

The command always begins with the CREST executable (`config.crest_executable`) followed by the temporary solute file. Flags are appended in the order defined by `MicrosolvatorConfig.build_flag_list()`.

## Custom Execution (`run_command`)

`Microsolvator.run` accepts a `run_command` callback if you need to redirect command execution (e.g., mocking in tests). The callback signature is:

```python
def run_command(command: Sequence[str], workdir: Path) -> subprocess.CompletedProcess[str]:
    ...
```

- `command` is the list described above.
- `workdir` is the directory containing the generated input files (solute.xyz, solvent.xyz, optional `.xcontrol`).
- The callback must return a `subprocess.CompletedProcess` with `stdout`, `stderr`, and `returncode` set.

When `run_command` is omitted, `subprocess.run` is used with `check=True`, `capture_output=True`, and `text=True`.

## Constraint Handling

If `constrain_solute=True` or `constrained_indices` are provided, an `.xcontrol` file is written in the working directory and `--nopreopt` is appended to the command (unless already enabled in the config).

## Working Directory Control

- `working_directory`: points microsolvator to an existing folder. Temporary files persist in this location.
- `keep_temps=True`: generates a unique temporary directory and keeps it after the run completes.
- Default (no `working_directory`, `keep_temps=False`): temporary files are discarded automatically.

The working directory always contains the files passed to CREST (`solute.xyz`, `solvent.xyz`, `.xcontrol` if any). Output files (`crest_best.xyz`, `full_ensemble.xyz`, `full_population.dat`) are read from the same directory.

## Result Introspection

`MicrosolvationResult` provides:

- `command`: list of shell tokens used to invoke CREST.
- `best_structure`: lowest-energy structure as an ASE `Atoms` instance, if produced.
- `ensemble`: list of `Atoms` for the entire ensemble (from `full_ensemble.xyz`).
- `population_path`: Path to `full_population.dat` when present.
- `stdout` / `stderr`: captured CREST output strings.

Use `result.ensure_outputs()` to assert that at least one structure was parsed.
