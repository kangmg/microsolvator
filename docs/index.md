# microsolvator

ASE-driven wrapper for CREST microsolvation and explicitly solvated reaction trajectory workflows.

---

## What it does

**`microsolvator`** runs CREST Quantum Cluster Growth (QCG) from Python, with ASE `Atoms` as inputs and outputs.

**`microsolvator.workflow`** takes a NEB/string-method trajectory and produces an explicitly solvated version ready for the next calculation.

```
gas-phase trajectory  →  CREST QCG  →  Packmol box  →  MD equilibration  →  solvated trajectory
```

---

## Quick start

### Microsolvation

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

config = MicrosolvatorConfig(
    nsolv=5,
    ensemble=True,
    implicit_model="alpb",
    implicit_solvent="h2o",
    threads=8,
)

result = Microsolvator.run(
    solute=read("ts.xyz"),
    solvent=read("water.xyz"),
    config=config,
    constrain_solute=True,
)

best = result.best_structure   # lowest-energy cluster (ensemble=True)
```

### Solvated trajectory

```python
from ase.io import read, write
from xtb.ase.calculator import XTB

from microsolvator import MicrosolvatorConfig
from microsolvator.workflow import (
    SolvatedTrajectoryBuilder,
    SolvationWorkflowConfig,
    PackmolConfig,
    EquilibrationConfig,
    RelaxationConfig,
)

images = read("neb_guess.traj", index=":")
water  = read("water.xyz")

config = SolvationWorkflowConfig(
    microsolv=MicrosolvatorConfig(nsolv=5, ensemble=True, threads=8),
    packmol=PackmolConfig(solvent_density=1.0, box_margin=12.0),
    equilibration=EquilibrationConfig(
        heating_schedule=[(100, 500), (200, 500), (300, 1000)],
        nvt_steps=5000,
    ),
    relaxation=RelaxationConfig(method="optimize", fmax=0.05),
)

result = SolvatedTrajectoryBuilder.build(
    reaction_images=images,
    solvent=water,
    config=config,
    calculator=XTB(method="GFN-FF"),
)

write("solvated_guess.traj", result.solvated_images)
```

---

## Package layout

| Package | Description |
|---------|-------------|
| `microsolvator` | CREST QCG wrapper — microsolvation of a single structure |
| `microsolvator.workflow` | Full pipeline — solvated reaction trajectory from a NEB guess |
