# Solvated Trajectory Workflow

`microsolvator.workflow` converts a gas-phase NEB / string-method trajectory into an explicitly solvated one.

## Pipeline

```
reaction_images[ts_index]
        │
        ▼
  1. CREST QCG          →  microsolvated cluster (solute + first shell)
        │
        ▼
  2. Packmol            →  cluster + bulk solvent in a periodic box
        │
        ▼
  3. MD equilibration   →  equilibrated box (solute frozen)
        │
        ▼
  4. Swap & relax       →  solvated image per reaction_image
```

Atom ordering is preserved throughout: `[solute | microsolv_shell | bulk_solvent]`.

---

## Quick example

```python
from ase.io import read, write
from xtb.ase.calculator import XTB
from microsolvator.workflow import solvate_trajectory

images = read("neb_guess.traj", index=":")
water  = read("water.xyz")

result = solvate_trajectory(
    images, water,
    calc=lambda: XTB(method="GFN-FF"),
    nsolv=5,
)

write("solvated_guess.traj", result.solvated_images)
```

`solvate_trajectory()` uses sensible defaults. For more control, pass a full `SolvationWorkflowConfig` via the `config=` parameter.

---

## Full example

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

# Load gas-phase trajectory
images = read("neb_guess.traj", index=":")
water  = read("water.xyz")

# Optional: set charge on input images
for img in images:
    img.info["charge"] = 0

config = SolvationWorkflowConfig(
    microsolv=MicrosolvatorConfig(
        nsolv=5,
        ensemble=True,             # best_structure is used
        implicit_model="alpb",
        implicit_solvent="h2o",
        threads=8,
    ),
    packmol=PackmolConfig(
        solvent_density=1.0,       # g/cm³
        box_margin=12.0,           # Å from cluster boundary to box wall
    ),
    equilibration=EquilibrationConfig(
        heating_schedule=[(100, 500), (200, 500), (300, 1000)],
        npt_steps=2000,            # optional NPT phase
        nvt_steps=5000,
        temperature=300.0,
    ),
    relaxation=RelaxationConfig(
        method="optimize",         # LBFGS
        fmax=0.05,
    ),
    ts_index=None,                 # auto: len(images) // 2
    working_directory="./solvation_run",
)

result = SolvatedTrajectoryBuilder.build(
    reaction_images=images,
    solvent=water,
    config=config,
    calc=lambda: XTB(method="GFN-FF"),
)

write("solvated_guess.traj", result.solvated_images)
```

---

## Calculator: instance vs factory

The `calc` parameter accepts either a calculator **instance** or a **factory** callable:

=== "Factory (recommended)"

    ```python
    # A fresh calculator is created for each MD/optimization step
    result = solvate_trajectory(
        images, water,
        calc=lambda: XTB(method="GFN-FF"),
        nsolv=5,
    )
    ```

=== "Instance"

    ```python
    # The instance is deep-copied internally for each step
    result = solvate_trajectory(
        images, water,
        calc=XTB(method="GFN-FF"),
        nsolv=5,
    )
    ```

!!! warning "Stateful calculators"
    Some calculators (GPAW, VASP, CP2K, …) carry internal state that doesn't survive
    `copy.deepcopy`. For these, **always use a factory**:

    ```python
    calc = lambda: GPAW(mode="fd", xc="PBE", txt=None)
    ```

---

## Practical examples

### S<sub>N</sub>2 reaction in water

```python
from ase.io import read, write
from ase.build import molecule
from xtb.ase.calculator import XTB
from microsolvator.workflow import solvate_trajectory

# Gas-phase NEB images for Cl⁻ + CH₃Br → ClCH₃ + Br⁻
images = read("sn2_neb.traj", index=":")

# Charge on the overall system
for img in images:
    img.info["charge"] = -1

water = molecule("H2O")

result = solvate_trajectory(
    images, water,
    calc=lambda: XTB(method="GFN2-xTB"),
    nsolv=8,              # 8 explicit water in first shell
    solvent_density=1.0,   # g/cm³ bulk water
    ts_index=4,            # known TS index
)

write("sn2_solvated.traj", result.solvated_images)
print(f"Total atoms per image: {result.n_total_atoms}")
print(f"  Solute: {result.n_solute_atoms}")
print(f"  Cluster (solute + shell): {result.n_cluster_atoms}")
```

### Diels-Alder in organic solvent

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

images = read("diels_alder_neb.traj", index=":")
toluene = read("toluene.xyz")

config = SolvationWorkflowConfig(
    microsolv=MicrosolvatorConfig(
        nsolv=4,
        ensemble=True,
        threads=4,
    ),
    packmol=PackmolConfig(
        solvent_density=0.87,      # toluene density
        box_margin=10.0,
    ),
    equilibration=EquilibrationConfig(
        heating_schedule=[(100, 300), (200, 300), (300, 500)],
        npt_steps=3000,            # NPT to relax box at target density
        nvt_steps=5000,
    ),
    relaxation=RelaxationConfig(
        method="fire",             # FIRE optimizer for stubborn geometries
        fmax=0.1,
        max_steps=300,
    ),
    working_directory="./da_solvation",
)

result = SolvatedTrajectoryBuilder.build(
    reaction_images=images,
    solvent=toluene,
    config=config,
    calc=lambda: XTB(method="GFN-FF"),
)

write("da_solvated.traj", result.solvated_images)
```

### Progress logging

```python
def log(step_name, info):
    print(f"[{step_name}] {info}")

result = solvate_trajectory(
    images, water,
    calc=lambda: XTB(method="GFN-FF"),
    nsolv=5,
    log_callback=log,
)
# Prints:
#   [microsolvation] {'ts_index': 4}
#   [packmol] {'n_cluster': 23}
#   [equilibration] {'n_fixed': 8}
#   [swap_relax] {'n_images': 9}
```

---

## Result fields

```python
result.solvated_images     # list[Atoms] — final solvated trajectory
result.microsolvated_ts    # Atoms — CREST QCG cluster (Step 1)
result.boxed_system        # Atoms — Packmol output (Step 2)
result.equilibrated_system # Atoms — MD template (Step 3)
result.n_solute_atoms      # int
result.n_cluster_atoms     # int — solute + microsolv shell
result.n_total_atoms       # int
result.ts_index            # int — resolved TS index
```

---

## Configuration reference

### `PackmolConfig`

| Field | Default | Description |
|-------|---------|-------------|
| `box_margin` | `10.0` | Å from cluster boundary to box wall |
| `solvent_density` | `None` | g/cm³ — auto-compute molecule count |
| `n_bulk_solvent` | `None` | explicit molecule count (alternative to density) |
| `min_distance` | `2.0` | Packmol tolerance (Å) |
| `packmol_executable` | `None` | path to packmol binary; None → PATH |
| `retry_margin_increment` | `2.0` | Å added to margin on each retry |
| `max_retries` | `3` | max Packmol retries before raising |

### `EquilibrationConfig`

| Field | Default | Description |
|-------|---------|-------------|
| `temperature` | `300.0` | K — NVT production and NPT temperature |
| `nvt_steps` | `5000` | NVT production steps |
| `npt_steps` | `0` | NPT steps; `0` disables NPT phase |
| `npt_pfactor` | `None` | barostat coupling (eV·fs²/Å³); None → 2.2 GPa estimate |
| `timestep` | `1.0` | fs |
| `friction` | `0.01` | Langevin friction (1/fs) |
| `heating_schedule` | `None` | `[(T_K, steps), …]` NVT ramp before production |

### `RelaxationConfig`

| Field | Default | Description |
|-------|---------|-------------|
| `method` | `"optimize"` | `"optimize"` (LBFGS), `"fire"`, or `"md"` (Langevin) |
| `fmax` | `0.05` | eV/Å convergence for optimize / fire |
| `max_steps` | `500` | max optimizer steps |
| `md_steps` | `200` | steps for `method="md"` |
| `md_temperature` | `300.0` | K for `method="md"` |
| `md_timestep` | `0.5` | fs for `method="md"` |
| `md_friction` | `0.01` | 1/fs for `method="md"` |

---

## Running individual steps

Each step is available as a standalone static method.

```python
from xtb.ase.calculator import XTB

calc = lambda: XTB(method="GFN-FF")

# Step 1 only
microsolv_result = SolvatedTrajectoryBuilder.microsolvate_ts(
    ts_image=images[3],
    solvent=water,
    config=config.microsolv,
    working_directory="./step1",
)
cluster = microsolv_result.best_structure or microsolv_result.final

# Step 2 only
boxed, n_bulk = SolvatedTrajectoryBuilder.pack_solvent_box(
    cluster=cluster,
    solvent=water,
    config=config.packmol,
)

# Step 3 only — freeze solute only; shell equilibrates with bulk
equilibrated = SolvatedTrajectoryBuilder.equilibrate(
    system=boxed,
    n_fixed=len(images[0]),       # solute atoms only
    calc=calc,
    config=config.equilibration,
    traj_path="equil.traj",
)

# Step 4 only
solvated_images = SolvatedTrajectoryBuilder.swap_and_relax(
    template=equilibrated,
    reaction_images=images,
    n_solute=len(images[0]),
    ts_index=3,
    calc=calc,
    config=config.relaxation,
)
```

---

## Charge and multiplicity

Set charge and multiplicity on the input images via `Atoms.info`. The builder propagates them through all steps and into the final images.

```python
for img in images:
    img.info["charge"] = -1
    img.info["multiplicity"] = 1
```

The `MicrosolvatorConfig.charge` / `.uhf` fields are overridden by these values when present.

---

## MD equilibration details

The equilibration runs in up to three phases:

1. **NVT heating** — stepwise ramp from `heating_schedule` using Langevin dynamics
2. **NPT** — optional; uses `ase.md.npt.NPT` (Nosé-Hoover barostat) — enable with `npt_steps > 0`
3. **NVT production** — Langevin dynamics at `temperature` for `nvt_steps` steps

Velocities are initialized from a Maxwell-Boltzmann distribution at the initial temperature before MD starts.

Only the solute atoms (indices `0 … n_solute-1`) are frozen with `ase.constraints.FixAtoms`.
The microsolvation shell is free to equilibrate together with the bulk solvent, forming a natural shell–bulk interface.

## Swap & alignment details

Starting from the TS image, the builder propagates **outward in both directions** (TS→product and TS→reactant).
Each image inherits the relaxed solvent from the previous (closer-to-TS) image as its starting template:

```
         ← TS-2 ← TS-1 ← [TS] → TS+1 → TS+2 →
```

For each non-TS image:

1. Swap the solute via Kabsch alignment onto the **previous image's** solute (not the original TS template)
2. Freeze only the solute (`0:n_solute`); microsolvation shell and bulk solvent are free to relax
3. Run constrained relaxation so the shell adapts to the new solute geometry

Each step gets a **fresh calculator** (via deep-copy or factory call) to prevent state leakage between images.

This chain propagation ensures:

- **Continuity** — adjacent images share similar solvent configurations
- **Shell adaptation** — the microsolvation shell re-optimises for each image's charge distribution and geometry

The TS image (`ts_index`) is returned as a direct copy of the equilibrated template — no swap, no drift.
