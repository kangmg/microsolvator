# microsolvator

ASE-driven wrapper for running CREST quantum cluster growth microsolvation workflows.

## Installation

```bash
pip install .
```

## Quick Start

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

solute = read("benzoic_acid.xyz")
solvent = read("water.xyz")

config = MicrosolvatorConfig(
    nsolv=3,
    implicit_model="alpb",
    implicit_solvent="h2o",
)

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrain_solute=True,
)

best = result.best_structure
ensemble = result.ensemble
```

## Implicit Solvent Check

```python
from microsolvator import supports_implicit_solvent

supports_implicit_solvent(method="gfn2", model="alpb", solvent="h2o")
```

## Execution Details

For an in-depth explanation of how commands are built, how to override execution with `run_command`, and how working directories are managed, see [`docs/usage.md`](docs/usage.md).
