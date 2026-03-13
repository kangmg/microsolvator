# Microsolvation

`Microsolvator.run` wraps CREST's Quantum Cluster Growth (QCG) mode. It takes ASE `Atoms` as inputs and returns parsed structures as `Atoms`.

## Basic usage

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

config = MicrosolvatorConfig(
    nsolv=5,
    method="gfn2",
    implicit_model="alpb",
    implicit_solvent="h2o",
    threads=8,
)

result = Microsolvator.run(
    solute=read("molecule.xyz"),
    solvent=read("water.xyz"),
    config=config,
)
```

## Result fields

| Field | File | When populated |
|-------|------|----------------|
| `result.best_structure` | `crest_best.xyz` | `ensemble=True` (default) |
| `result.final` | `grow/cluster.xyz` | Always (grow phase) |
| `result.ensemble` | `full_ensemble.xyz` | `ensemble=True` |
| `result.traj` | `grow/qcg_grow.xyz` | Always (grow phase) |

**Which one to use:**

- `ensemble=True` → use `result.best_structure` (conformer-searched, lower energy)
- `ensemble=False` → use `result.final` (grow result only)

## Constraining the solute

```python
result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrain_solute=True,        # freeze all solute atoms
    # constrained_indices=[1, 2],  # or specify 1-based indices
)
```

A `.xcontrol` file is written automatically; `--nopreopt` is appended to the CREST command.

## Keeping files / specifying working directory

```python
from pathlib import Path

# Persist files in a specific directory
result = Microsolvator.run(
    solute=solute, solvent=solvent, config=config,
    working_directory=Path("./my_run"),
)

# Keep auto-generated temp directory
result = Microsolvator.run(
    solute=solute, solvent=solvent, config=config,
    keep_temps=True,
)

print(result.working_directory)  # Path to temp dir
```

## Dry run

Inspect the CREST command without executing it.

```python
result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    working_directory=Path("./scratch"),
    prepare_only=True,
)

print(result.shell_command)
# /path/to/crest solute.xyz --qcg solvent.xyz --nsolv 5 ...
```

## Implicit solvent support

```python
from microsolvator import supports_implicit_solvent, list_supported_implicit_solvents

supports_implicit_solvent(method="gfn2", model="alpb", solvent="h2o")  # True

# List all solvents for a given model
list_supported_implicit_solvents(method="gfn2", model="alpb")
```

## Logging

By default, CREST output is streamed to `crest_run.log` in the working directory.

```python
result = Microsolvator.run(
    solute=solute, solvent=solvent, config=config,
    log_file="my_run.log",  # custom path
    # log_file=None,         # disable logging
)
```

## Custom executor

```python
import subprocess
from pathlib import Path

def my_runner(command, workdir, env=None):
    # e.g., submit to a cluster scheduler
    print("Would run:", " ".join(command))
    return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

result = Microsolvator.run(
    solute=solute, solvent=solvent, config=config,
    run_command=my_runner,
)
```
