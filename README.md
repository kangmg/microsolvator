# microsolvator

ASE-driven wrapper for running CREST quantum cluster growth microsolvation workflows.

## Installation

```bash
pip install git+https://github.com/kangmg/microsolvator.git
```

## Quick Start

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
    crest_executable="/abs/path/crest",
    xtb_executable="/abs/path/xtb",
    threads=12,
)

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrain_solute=True,
    log_file="crest_run.log",
)

best = result.best_structure
ensemble = result.ensemble
final_cluster = result.final
trajectory = result.traj

# Dry run to inspect the command and generated inputs
dry_run = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    working_directory=Path("./scratch"),
    prepare_only=True,
)

print(dry_run.shell_command)
```

`threads` maps directly to CREST's `--T` (number of threads); temperature is now passed via `--temp`.

## Implicit Solvent Check

```python
from microsolvator import supports_implicit_solvent

supports_implicit_solvent(method="gfn2", model="alpb", solvent="h2o")
```

## Execution Details

For an in-depth explanation of how commands are built, how to override execution with `run_command`, and how working directories are managed, see [`docs/usage.md`](docs/usage.md).

### Binary Resolution

`Microsolvator` resolves CREST and xTB executables using, in order, explicit config paths, environment variables (`CREST_BIN`, `XTB_BIN`), binaries bundled in `microsolvator/_bin/` (installable via `install_crest`/`install_xtb`), and finally the system `PATH`.

Resolved paths are injected into the command line as absolute paths (including an explicit `--xnam /abs/path/to/xtb` flag) when available.

### Logging

During execution, CREST stdout/stderr are streamed in real time to `crest_run.log` (relative to the working directory by default). Override with `log_file="my_crest.log"` or disable by passing `log_file=None`.
