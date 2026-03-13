# microsolvator.workflow

## `SolvatedTrajectoryBuilder`

```python
from microsolvator.workflow import SolvatedTrajectoryBuilder
```

### `build` — full pipeline

```python
SolvatedTrajectoryBuilder.build(
    *,
    reaction_images: list[Atoms],
    solvent: Atoms,
    config: SolvationWorkflowConfig,
    calculator,
    log_callback: Callable[[str, Any], None] | None = None,
) -> SolvatedTrajectoryResult
```

| Parameter | Description |
|-----------|-------------|
| `reaction_images` | NEB / string-method images |
| `solvent` | One solvent molecule |
| `config` | `SolvationWorkflowConfig` |
| `calculator` | Any ASE-compatible calculator |
| `log_callback` | Called at each step boundary: `callback(step_name, info)` |

### Individual steps

```python
# Step 1
SolvatedTrajectoryBuilder.microsolvate_ts(ts_image, solvent, config, **kwargs)
    -> MicrosolvationResult

# Step 2
SolvatedTrajectoryBuilder.pack_solvent_box(cluster, solvent, config, workdir=None)
    -> tuple[Atoms, int]  # (boxed_system, n_bulk_solvent)

# Step 3
SolvatedTrajectoryBuilder.equilibrate(system, n_fixed, calculator, config,
                                       traj_path=None, log_callback=None)
    -> Atoms

# Step 4
SolvatedTrajectoryBuilder.swap_and_relax(template, reaction_images, n_solute,
                                          ts_index, calculator, config)
    -> list[Atoms]
```

---

## `SolvationWorkflowConfig`

```python
from microsolvator.workflow import SolvationWorkflowConfig
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `microsolv` | `MicrosolvatorConfig` | *(required)* | CREST QCG config for Step 1 |
| `packmol` | `PackmolConfig` | `PackmolConfig()` | Packmol config for Step 2 |
| `equilibration` | `EquilibrationConfig` | `EquilibrationConfig()` | MD config for Step 3 |
| `relaxation` | `RelaxationConfig` | `RelaxationConfig()` | Relax config for Step 4 |
| `ts_index` | `int \| None` | `None` | TS image index; `None` → `len(images) // 2` |
| `working_directory` | `Path \| None` | `None` | Intermediate files directory |
| `keep_temps` | `bool` | `False` | Keep auto-generated temp dir |

---

## `PackmolConfig`

```python
from microsolvator.workflow import PackmolConfig
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `box_margin` | `float` | `10.0` | Å from cluster to box wall |
| `solvent_density` | `float \| None` | `None` | g/cm³ — auto-compute molecule count |
| `n_bulk_solvent` | `int \| None` | `None` | Explicit molecule count |
| `min_distance` | `float` | `2.0` | Packmol tolerance (Å) |
| `packmol_executable` | `str \| None` | `None` | Packmol binary path; `None` → PATH |
| `retry_margin_increment` | `float` | `2.0` | Å added on each retry |
| `max_retries` | `int` | `3` | Max retries before raising |

Provide exactly one of `solvent_density` or `n_bulk_solvent`.

---

## `EquilibrationConfig`

```python
from microsolvator.workflow import EquilibrationConfig
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `temperature` | `float` | `300.0` | K — NVT/NPT target temperature |
| `nvt_steps` | `int` | `5000` | NVT production steps (Langevin) |
| `npt_steps` | `int` | `0` | NPT steps; `0` skips NPT phase |
| `npt_pfactor` | `float \| None` | `None` | NPT barostat coupling (eV·fs²/Å³); `None` → 2.2 GPa estimate |
| `timestep` | `float` | `1.0` | fs |
| `friction` | `float` | `0.01` | Langevin friction (1/fs) |
| `heating_schedule` | `list[tuple[float, int]] \| None` | `None` | `[(T_K, steps), …]` |

MD phases: NVT heating → NPT (optional) → NVT production.

---

## `RelaxationConfig`

```python
from microsolvator.workflow import RelaxationConfig
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `method` | `str` | `"optimize"` | `"optimize"` (LBFGS), `"fire"`, or `"md"` (Langevin) |
| `fmax` | `float` | `0.05` | eV/Å convergence threshold |
| `max_steps` | `int` | `500` | Max optimizer steps |
| `md_steps` | `int` | `200` | Steps for `method="md"` |
| `md_temperature` | `float` | `300.0` | K for `method="md"` |
| `md_timestep` | `float` | `0.5` | fs for `method="md"` |
| `md_friction` | `float` | `0.01` | 1/fs for `method="md"` |

---

## `SolvatedTrajectoryResult`

```python
from microsolvator.workflow import SolvatedTrajectoryResult
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `solvated_images` | `list[Atoms]` | Final solvated trajectory |
| `microsolvated_ts` | `Atoms` | CREST cluster (Step 1) |
| `boxed_system` | `Atoms` | Packmol output with cell (Step 2) |
| `equilibrated_system` | `Atoms` | MD-equilibrated template (Step 3) |
| `n_solute_atoms` | `int` | Atoms in one NEB image |
| `n_cluster_atoms` | `int` | Solute + microsolvation shell atoms |
| `n_total_atoms` | `int` | Total atoms in solvated box |
| `ts_index` | `int` | Resolved TS index |
| `config` | `SolvationWorkflowConfig` | Config used |
| `working_directory` | `Path` | Where intermediate files were written |

Atom ordering in all images: `[solute | microsolv_shell | bulk_solvent]`.

---

## Low-level utilities

```python
from microsolvator.workflow.swap import kabsch_align, swap_solute, relax_interface
from microsolvator.workflow.packmol import PackmolSolvator
from microsolvator.workflow.equilibration import equilibrate
from microsolvator.workflow.utils import (
    estimate_box_size,
    count_solvent_molecules,
    validate_packmol_output,
)
```
