# Solvated Reaction Trajectory Builder — Workflow Design

## Problem Statement

NEB / string method 계산을 위한 초기 guess trajectory가 주어졌을 때,
explicit solvation 환경을 구축하여 **용매화된 반응 경로 trajectory guess**를 자동으로 생성한다.

## Input / Output

### Input
- **reaction_images**: `list[Atoms]` — NEB/string method의 초기 guess trajectory (reactant → TS → product)
- **solvent**: `Atoms` — 용매 분자 1개 (e.g., water)
- **config**: `SolvationWorkflowConfig` — 전체 워크플로우 설정

### Output
- **SolvatedTrajectoryResult**: 용매화된 전체 reaction trajectory (`list[Atoms]`)
  - 각 이미지가 동일한 solvation box 내에 위치
  - 반응물/생성물 영역만 달라지고, 용매 환경은 공유

---

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    Solvated Trajectory Builder                       │
│                                                                     │
│  ┌──────────┐    ┌──────────────┐    ┌──────────┐    ┌───────────┐ │
│  │ 1. CREST │    │ 2. Packmol   │    │ 3. MD    │    │ 4. Swap   │ │
│  │ Micro-   │───▶│ Box          │───▶│ Equili-  │───▶│ Endpoints │ │
│  │ solvation│    │ Solvation    │    │ bration  │    │ & Relax   │ │
│  └──────────┘    └──────────────┘    └──────────┘    └───────────┘ │
│       │                                                    │        │
│   TS image                                          Solvated traj  │
│  (middle)                                           (all images)    │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Step 1: Microsolvation of TS (CREST QCG)

### 목적
TS guess 구조 주변에 first solvation shell을 정밀하게 배치한다.
CREST의 Quantum Cluster Growth (QCG) 방법으로 energetically favorable한 위치에 용매를 배치.

### CREST QCG 결과 해석

CREST QCG는 두 단계로 동작한다.

| 단계 | 파일 | `MicrosolvationResult` 필드 | 의미 |
|------|------|--------------------------|------|
| Grow | `grow/cluster.xyz` | `result.final` | 용매 배치 완료된 클러스터 |
| Ensemble | `crest_best.xyz` | `result.best_structure` | 앙상블 탐색에서 에너지 최저 conformer |
| Ensemble | `full_ensemble.xyz` | `result.ensemble` | 전체 conformer 목록 |
| Grow | `grow/qcg_grow.xyz` | `result.traj` | 용매 추가 과정 trajectory |

**결과 선택 로직**:
- `ensemble=True` (기본값, `MicrosolvatorConfig.ensemble`): Grow + Ensemble 수행
  → `result.best_structure` 사용 (더 정제된 결과)
- `ensemble=False`: Grow만 수행
  → `result.final` 사용 (`best_structure`는 None)

**올바른 fallback 순서**: `result.best_structure or result.final`
(반대 순서는 ensemble=True일 때 더 정제된 결과를 버리는 것이므로 잘못됨)

### 상세

```python
# pseudo-code
ts_image = reaction_images[ts_index]  # ts_index: 명시 또는 len(images) // 2

microsolv_result = Microsolvator.run(
    solute=ts_image,
    solvent=solvent,
    config=microsolv_config,      # MicrosolvatorConfig.ensemble 에 따라 결과가 달라짐
    constrain_solute=True,
)

# ensemble=True  → best_structure (앙상블 탐색 결과)
# ensemble=False → final          (grow 결과)
microsolvated_ts = microsolv_result.best_structure or microsolv_result.final
if microsolvated_ts is None:
    raise ValueError("Microsolvation produced no valid cluster structure")

n_solute = len(ts_image)
n_cluster = len(microsolvated_ts)   # solute + microsolv shell 원자 수
```

### 설계 포인트
- `ts_index` 사용자 지정 가능. 기본값은 `len(images) // 2`
- `MicrosolvatorConfig.ensemble`이 결과 선택 로직을 결정하므로 config와 결과 선택 로직을 연동
- Charge/multiplicity: `Atoms.info["charge"]`와 `Atoms.info["multiplicity"]` 컨벤션으로 전달 (Step 3, 4 동일)

---

## Step 2: Packmol Box Solvation

### 목적
Microsolvated cluster를 box 중심에 놓고, 주변을 bulk 용매로 채운 시뮬레이션 박스를 생성한다.

### 상세

```python
# pseudo-code
box_system = PackmolSolvator.run(
    solute_cluster=microsolvated_ts,
    solvent=solvent,
    config=packmol_config,
)
# box_system.cell → simulation box vectors
# atoms[0:n_cluster] = microsolvated cluster (순서 보존)
# atoms[n_cluster:]  = bulk solvent molecules
```

**전처리: cluster를 box center로 이동**

Packmol `fixed` 키워드는 입력 좌표 기준으로 위치를 고정한다.
Cluster가 임의의 좌표에 있을 수 있으므로, packmol input 생성 전에
**cluster를 box center로 이동**해야 box 안에 올바르게 배치된다.

```python
# Cluster를 box center로 translation (회전 없음)
box_center = np.array([box_size / 2] * 3)
cluster_com = microsolvated_ts.get_center_of_mass()
translation = box_center - cluster_com
centered_cluster = microsolvated_ts.copy()
centered_cluster.translate(translation)

# Packmol input: cluster는 fixed (0. 0. 0. 0. 0. 0. → 이미 center로 이동했으므로)
# structure cluster.xyz
#   number 1
#   fixed 0. 0. 0. 0. 0. 0.   ← 좌표 그대로 (이미 center에 위치)
#   center                     ← packmol center 키워드 필요 없음
# end structure
```

**Packmol output 후 원자 순서 검증**

```python
n_solvent_atoms = len(solvent)
expected_total = n_cluster + n_bulk_solvent * n_solvent_atoms
assert len(box_system) == expected_total, (
    f"Packmol output atom count mismatch: "
    f"expected {expected_total}, got {len(box_system)}"
)
# 첫 n_cluster 원자가 cluster와 동일한지 원소 검증
assert list(box_system.symbols[:n_cluster]) == list(centered_cluster.symbols), \
    "Atom order in Packmol output does not match cluster"
```

### 설계 포인트
- 용매 분자 수: `PackmolConfig.solvent_density` (g/cm³) 기반 자동 계산 또는 `n_bulk_solvent` 직접 지정
- Packmol 실패 시 (`n_bulk_solvent`가 너무 많거나 box가 너무 작으면) `box_margin`을 키워서 자동 retry
- 출력 Atoms의 원자 순서: `[solute | microsolv_shell | bulk_solvent]` — 이 순서를 검증 후 보장

---

## Step 3: MD Equilibration (Constrained, 2-Phase)

### 목적
Bulk 용매를 평형화하되, microsolvated cluster (solute + first shell)는 고정한다.

### 2-Phase 평형화

1단계 NVT만으로는 Packmol이 생성한 초기 구조의 density artifact가 남을 수 있으므로,
**NVT heating → NPT density relaxation → NVT production** 3단계를 기본으로 한다.

```python
# pseudo-code
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin
from ase.md.npt import NPT

system = boxed_system.copy()
# Cluster (solute + microsolv shell) 고정
cluster_constraint = FixAtoms(indices=list(range(n_cluster)))
system.set_constraint(cluster_constraint)
system.calc = calculator  # 사용자 주입 ASE calculator

# Phase 1: NVT heating
if config.equilibration.heating_schedule:
    for temp_K, steps in config.equilibration.heating_schedule:
        dyn = Langevin(system, timestep=config.equilibration.timestep * units.fs,
                       temperature_K=temp_K, friction=config.equilibration.friction)
        dyn.run(steps)

# Phase 2: NPT density relaxation (선택)
if config.equilibration.npt_steps > 0:
    dyn = NPT(system, timestep=config.equilibration.timestep * units.fs,
              temperature_K=config.equilibration.temperature,
              externalstress=0.0, ttime=25 * units.fs, pfactor=...)
    dyn.run(config.equilibration.npt_steps)

# Phase 3: NVT production
dyn = Langevin(system, timestep=config.equilibration.timestep * units.fs,
               temperature_K=config.equilibration.temperature,
               friction=config.equilibration.friction)
dyn.run(config.equilibration.nvt_steps)

equilibrated_system = system.copy()
```

### Charge/multiplicity 전달

```python
charge = boxed_system.info.get("charge", 0)
mult   = boxed_system.info.get("multiplicity", 1)
# calculator에 charge/multiplicity 반영은 calculator 종류에 따라 다름
# → EquilibrationConfig.charge_key, multiplicity_key로 설정 가능하게
```

### 설계 포인트
- PBC: cell이 설정되어 있으면 `system.pbc = True` 자동 설정
- Trajectory 저장: `equilibration_traj.traj` (ASE Trajectory)
- Calculator는 사용자가 주입 (calculator-agnostic)

---

## Step 4: Endpoint Swap & Relaxation

### 목적
평형화된 용매 환경을 유지하면서, 반응물/생성물 구조로 solute만 교체하고 interface를 relaxation한다.

### 4a. TS 이미지 처리

**TS 이미지 (ts_index)는 swap을 건너뛴다.**
equilibrated_system의 solute가 이미 TS이므로, 다시 swap + relax하면 수치적 drift만 발생한다.

```python
solvated_images = []
for i, image in enumerate(reaction_images):
    if i == ts_index:
        # TS는 equilibrated system을 그대로 사용
        solvated_images.append(equilibrated_system.copy())
        continue

    swapped = swap_solute(equilibrated_system, image, n_solute)
    swapped = relax_interface(swapped, n_solute, calculator, config.relaxation)
    solvated_images.append(swapped)
```

### 4b. Geometry Swap — Kabsch Alignment

COM alignment만으로는 orientation을 맞출 수 없다.
**Kabsch alignment** (RMSD 최소화하는 rotation + translation)를 기본으로 사용한다.
ASE 내장 `minimize_rotation_and_translation`을 활용하면 scipy 없이 구현 가능.

```python
from ase.geometry import get_distances

def swap_solute(solvated_system: Atoms, new_solute: Atoms, n_solute: int) -> Atoms:
    """Replace solute atoms with Kabsch-aligned new geometry."""
    new_system = solvated_system.copy()

    ref_solute = solvated_system[:n_solute].copy()   # 기준: equilibrated TS solute
    aligned_solute = new_solute.copy()

    # ASE Kabsch-like alignment: ref에 aligned를 맞춤
    minimize_rotation_and_translation(ref_solute, aligned_solute)

    new_system.positions[:n_solute] = aligned_solute.positions
    return new_system
```

> `minimize_rotation_and_translation(ref, atoms)` — `atoms`를 `ref`에 맞게 in-place 회전+평행이동.
> 반응 경로에서 원소 순서가 보존되어야 하며, 원소 불일치 시 ValueError를 명시적으로 발생시킨다.

### 4c. Constrained Relaxation

```python
def relax_interface(system: Atoms, n_solute: int, calculator, config: RelaxationConfig) -> Atoms:
    system = system.copy()
    # Solute 고정, 용매(microsolv + bulk)만 relax
    system.set_constraint(FixAtoms(indices=list(range(n_solute))))
    system.calc = calculator

    if config.method == "optimize":
        opt = LBFGS(system, logfile=None)
        opt.run(fmax=config.fmax, steps=config.max_steps)
    elif config.method == "fire":
        opt = FIRE(system, logfile=None)
        opt.run(fmax=config.fmax, steps=config.max_steps)
    elif config.method == "md":
        dyn = Langevin(system, timestep=config.md_timestep * units.fs,
                       temperature_K=config.md_temperature, friction=config.md_friction)
        dyn.run(config.md_steps)

    return system
```

---

## Configuration Design (분리된 구조)

단일 monolithic config 대신, step별 config를 분리하여 **개별 step 독립 호출 시 불필요한 필드가 없도록** 설계한다.

```python
@dataclass
class PackmolConfig:
    box_margin: float = 10.0               # Å, cluster boundary로부터의 여유
    solvent_density: Optional[float] = None # g/cm³, None이면 n_bulk_solvent 사용
    n_bulk_solvent: Optional[int] = None    # 직접 지정 (density와 둘 중 하나)
    min_distance: float = 2.0              # Å, 분자 간 최소 거리
    packmol_executable: Optional[str] = None  # None이면 PATH에서 검색
    retry_margin_increment: float = 2.0    # Å, retry 시 box_margin 증가량
    max_retries: int = 3


@dataclass
class EquilibrationConfig:
    temperature: float = 300.0             # K, NVT production temperature
    nvt_steps: int = 5000                  # NVT production steps
    npt_steps: int = 2000                  # NPT density relaxation steps (0이면 skip)
    timestep: float = 1.0                  # fs
    friction: float = 0.01                 # 1/fs (Langevin)
    heating_schedule: Optional[list[tuple[float, int]]] = None  # [(T_K, steps), ...]
    # e.g. [(100, 500), (200, 500), (300, 1000)] → 단계적 heating


@dataclass
class RelaxationConfig:
    method: str = "optimize"               # "optimize" | "fire" | "md"
    fmax: float = 0.05                     # eV/Å (optimize, fire)
    max_steps: int = 500
    md_steps: int = 200                    # method="md" 전용
    md_temperature: float = 300.0          # K, method="md" 전용
    md_timestep: float = 0.5              # fs, method="md" 전용
    md_friction: float = 0.01             # 1/fs, method="md" 전용


@dataclass
class SolvationWorkflowConfig:
    """전체 워크플로우 설정 — step별 config를 조합"""

    # Step 1: Microsolvation
    microsolv: MicrosolvatorConfig         # MicrosolvatorConfig.ensemble 이 결과 선택을 결정

    # Step 2: Packmol
    packmol: PackmolConfig = field(default_factory=PackmolConfig)

    # Step 3: MD Equilibration
    equilibration: EquilibrationConfig = field(default_factory=EquilibrationConfig)

    # Step 4: Relaxation
    relaxation: RelaxationConfig = field(default_factory=RelaxationConfig)

    # General
    ts_index: Optional[int] = None         # None이면 자동 (len(images) // 2)
    working_directory: Optional[Path] = None
    keep_temps: bool = False
```

---

## Charge/Multiplicity 전달 체계

`Atoms.info` 딕셔너리를 통해 일관된 방식으로 전달한다.

```python
# 사용자가 입력 시 설정
ts_image.info["charge"] = -1
ts_image.info["multiplicity"] = 1

# Step 1: MicrosolvatorConfig.charge / .uhf 로 변환 (워크플로우 내부에서 자동)
# Step 2: Packmol은 charge 무관 (topology만 다루므로)
# Step 3, 4: calculator에 charge 전달 방식은 calculator-agnostic
#   → charge/multiplicity를 calculator 설정 시 사용자가 직접 반영하거나,
#     WorkflowConfig에 charge: int = 0, multiplicity: int = 1 필드를 두고
#     builder가 calculator wrapper를 통해 주입

# builder 내부에서 charge/mult를 Atoms.info에서 읽어 microsolv_config에 반영하는 로직
charge = reaction_images[ts_index].info.get("charge", config.microsolv.charge)
# MicrosolvatorConfig.charge는 여전히 직접 지정도 가능 (기존 API 호환)
```

**원칙**:
- `MicrosolvatorConfig.charge` / `.uhf` — CREST 전달용 (기존 필드 유지)
- `Atoms.info["charge"]` / `Atoms.info["multiplicity"]` — 구조 간 전달 컨벤션
- Step 3/4에서의 calculator charge 반영은 **사용자 책임** (calculator-agnostic 원칙 유지)
- builder는 `Atoms.info`의 값을 `MicrosolvatorConfig`에 override하는 convenience logic 제공

---

## Public API Design

```python
class SolvatedTrajectoryBuilder:
    """Build explicitly solvated reaction trajectories."""

    @classmethod
    def build(
        cls,
        *,
        reaction_images: list[Atoms],
        solvent: Atoms,
        config: SolvationWorkflowConfig,
        calculator,                        # ASE calculator (pluggable)
        log_callback: Optional[Callable[[str, Any], None]] = None,
    ) -> SolvatedTrajectoryResult:
        """Full pipeline: microsolvate → pack → equilibrate → swap & relax."""
        ...

    # 개별 step 독립 호출 — step별 config만 전달
    @staticmethod
    def microsolvate_ts(
        ts_image: Atoms,
        solvent: Atoms,
        config: MicrosolvatorConfig,
        **kwargs,
    ) -> MicrosolvationResult:
        """Step 1 only."""
        ...

    @staticmethod
    def pack_solvent_box(
        cluster: Atoms,
        solvent: Atoms,
        config: PackmolConfig,
    ) -> Atoms:
        """Step 2 only."""
        ...

    @staticmethod
    def equilibrate(
        system: Atoms,
        n_fixed: int,
        calculator,
        config: EquilibrationConfig,
    ) -> Atoms:
        """Step 3 only."""
        ...

    @staticmethod
    def swap_and_relax(
        template: Atoms,
        reaction_images: list[Atoms],
        n_solute: int,
        ts_index: int,
        calculator,
        config: RelaxationConfig,
    ) -> list[Atoms]:
        """Step 4 only. ts_index 이미지는 template을 그대로 반환."""
        ...


@dataclass
class SolvatedTrajectoryResult:
    """워크플로우 전체 결과"""
    solvated_images: list[Atoms]           # 최종 용매화된 trajectory
    microsolvated_ts: Atoms                # Step 1 결과
    boxed_system: Atoms                    # Step 2 결과
    equilibrated_system: Atoms             # Step 3 결과
    n_solute_atoms: int
    n_cluster_atoms: int                   # solute + microsolv shell
    n_total_atoms: int
    ts_index: int                          # 실제 사용된 ts_index
    config: SolvationWorkflowConfig
    working_directory: Path
```

---

## Package Structure

```
microsolvator/
├── __init__.py                    # 기존 + 신규 export
├── config.py                     # 기존 MicrosolvatorConfig
├── command.py                    # 기존 CREST command builder
├── runner.py                     # 기존 Microsolvator
├── results.py                    # 기존 MicrosolvationResult
├── support.py                    # 기존 implicit solvent support
├── install.py                    # 기존 binary installer
│
├── workflow/                     # ★ 신규 서브패키지
│   ├── __init__.py
│   ├── config.py                 # PackmolConfig, EquilibrationConfig, RelaxationConfig, SolvationWorkflowConfig
│   ├── builder.py                # SolvatedTrajectoryBuilder (main orchestrator)
│   ├── results.py                # SolvatedTrajectoryResult
│   ├── packmol.py                # Packmol wrapper (input gen, execution, parsing, validation)
│   ├── equilibration.py          # MD equilibration (NVT heating → NPT → NVT)
│   ├── swap.py                   # swap_solute (Kabsch alignment), relax_interface
│   └── utils.py                  # box size/density calculation, atom order validation
│
tests/
├── test_microsolvator.py         # 기존
├── test_workflow/                # ★ 신규
│   ├── test_packmol.py
│   ├── test_equilibration.py
│   ├── test_swap.py
│   └── test_builder.py           # integration test
```

---

## Dependency 정리

| Component | 의존성 | 비고 |
|-----------|--------|------|
| Step 1 (Microsolvation) | `crest`, `xtb` | 기존 microsolvator |
| Step 2 (Packmol) | `packmol` | 외부 바이너리 |
| Step 3 (MD) | ASE | calculator는 사용자 주입 |
| Step 4 (Swap) | ASE (`minimize_rotation_and_translation`) | scipy 불필요 |
| Step 4 (Relax) | ASE | calculator는 사용자 주입 |

`pyproject.toml` 추가 의존성: 없음 (scipy 불필요, ASE 내장 활용)

---

## Usage Example

```python
from ase.io import read, write
from microsolvator import MicrosolvatorConfig
from microsolvator.workflow import (
    SolvatedTrajectoryBuilder,
    SolvationWorkflowConfig,
    PackmolConfig,
    EquilibrationConfig,
    RelaxationConfig,
)

# 1. 반응 경로 로드
reaction_images = read("neb_initial_guess.traj", index=":")
solvent = read("water.xyz")

# Charge 설정 (필요 시)
for img in reaction_images:
    img.info["charge"] = 0

# 2. Calculator 준비 (사용자 선택)
from xtb.ase.calculator import XTB
calculator = XTB(method="GFN-FF")

# 3. 워크플로우 설정
config = SolvationWorkflowConfig(
    microsolv=MicrosolvatorConfig(
        nsolv=5,
        ensemble=True,             # True → best_structure, False → final
        implicit_model="alpb",
        implicit_solvent="h2o",
        threads=8,
    ),
    packmol=PackmolConfig(
        box_margin=12.0,
        solvent_density=1.0,
    ),
    equilibration=EquilibrationConfig(
        heating_schedule=[(100, 500), (200, 500), (300, 1000)],
        npt_steps=2000,
        nvt_steps=5000,
        temperature=300.0,
    ),
    relaxation=RelaxationConfig(
        method="optimize",
        fmax=0.05,
    ),
)

# 4. 실행
result = SolvatedTrajectoryBuilder.build(
    reaction_images=reaction_images,
    solvent=solvent,
    config=config,
    calculator=calculator,
)

# 5. 결과
write("solvated_neb_guess.traj", result.solvated_images)
```

---

## Edge Cases & 고려사항

1. **원자 순서 일관성**: 모든 step에서 `[solute | microsolv | bulk]` 순서 유지. Packmol output 후 반드시 검증.
2. **Charged systems**: `Atoms.info["charge"]` 컨벤션으로 전달. CREST로는 `MicrosolvatorConfig.charge` 자동 반영.
3. **Large conformational change**: reactant↔product 구조 차이가 크면 Kabsch alignment가 실패할 수 있음 → RMSD warning 제공.
4. **Packmol failure**: box가 너무 작으면 `box_margin`을 `retry_margin_increment`씩 늘려서 최대 `max_retries`회 재시도.
5. **PBC consistency**: equilibrated_system의 cell 정보가 모든 solvated_images에 복사되어야 함.
6. **TS index out of bounds**: `ts_index >= len(reaction_images)` 시 ValueError.
7. **Single-image trajectory**: NEB가 아닌 단일 구조 용매화에도 대응 가능하도록 `len(images) == 1` 처리.
