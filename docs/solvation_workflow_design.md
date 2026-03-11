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
- **SolvatedTrajectory**: 용매화된 전체 reaction trajectory (`list[Atoms]`)
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

### 상세
- **입력**: `reaction_images[ts_index]` (기본: 중간 이미지), `solvent`
- **처리**: 기존 `microsolvator.Microsolvator.run()` 호출
  - `constrain_solute=True`: solute geometry 보존
  - `nsolv`: first shell 용매 분자 수 (사용자 지정 또는 자동 추정)
- **출력**: `microsolvated_cluster: Atoms` (solute + n개 용매)
- **기록**: `n_solute_atoms`, `n_microsolv_atoms` (constraint에 사용)

```python
# pseudo-code
ts_image = reaction_images[len(reaction_images) // 2]
microsolv_result = Microsolvator.run(
    solute=ts_image,
    solvent=solvent,
    config=microsolv_config,
    constrain_solute=True,
)
microsolvated_ts = microsolv_result.best_structure  # solute + first shell
n_solute = len(ts_image)
n_cluster = len(microsolvated_ts)  # solute + microsolv shell
```

### 설계 포인트
- `ts_index` 를 사용자가 지정 가능. 기본값은 `len(images) // 2`
- `nsolv` 자동 추정 옵션: 용질의 solvent-accessible surface area 기반 (optional, v2)
- microsolvation 결과에서 `best_structure` 사용 (lowest energy cluster)

---

## Step 2: Packmol Box Solvation

### 목적
Microsolvated cluster를 중심에 놓고, 주변을 bulk 용매로 채운 시뮬레이션 박스를 생성한다.

### 상세
- **입력**: `microsolvated_ts`, `solvent`, box parameters
- **처리**:
  1. Microsolvated cluster의 bounding box 계산
  2. 여유 margin을 더해 시뮬레이션 box 크기 결정
  3. Packmol input 생성: cluster는 fixed, 나머지 공간에 solvent 채움
  4. Packmol 실행
- **출력**: `boxed_system: Atoms` (cluster + bulk solvent, with cell info)

```python
# pseudo-code
box_system = PackmolSolvator.run(
    solute_cluster=microsolvated_ts,
    solvent=solvent,
    box_margin=10.0,          # Å margin from cluster boundary
    solvent_density=1.0,      # g/cm³ (or n_solvent directly)
    min_distance=2.0,         # Å, min distance between molecules
    packmol_executable="packmol",
)
# box_system.cell → simulation box vectors
# atoms[0:n_cluster] = microsolvated cluster (순서 보존)
# atoms[n_cluster:] = bulk solvent
```

### 설계 포인트
- Packmol input file을 자동 생성 (`structure`, `inside box`, `fixed` 키워드)
- Cluster를 box 중심에 고정 배치 (packmol의 `fixed` 키워드)
- 용매 분자 수: density 기반 자동 계산 또는 직접 지정
- Packmol executable 경로: config 또는 PATH에서 resolve
- 출력 Atoms의 원자 순서: `[solute | microsolv_shell | bulk_solvent]` — 이 순서 보장이 핵심

---

## Step 3: MD Equilibration (Constrained)

### 목적
Bulk 용매를 평형화하되, microsolvated cluster (solute + first shell)는 고정한다.

### 상세
- **입력**: `boxed_system`, constraint indices (`0 ~ n_cluster-1`)
- **처리**:
  1. ASE `FixAtoms` constraint로 cluster 원자들 고정
  2. Calculator 설정 (GFN-FF 또는 사용자 지정 ASE calculator)
  3. MD 실행 (NVT 또는 NPT)
     - 초기 heating → equilibration
  4. 최종 프레임 추출
- **출력**: `equilibrated_system: Atoms`

```python
# pseudo-code
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin

system = boxed_system.copy()
system.set_constraint(FixAtoms(indices=list(range(n_cluster))))

# Calculator 설정 (pluggable)
system.calc = calculator  # e.g., XTB, GFN-FF, or any ASE calculator

# MD 실행
dyn = Langevin(system, timestep=1.0*units.fs, temperature_K=300, friction=0.01)
dyn.run(steps=md_steps)

equilibrated_system = system.copy()
```

### 설계 포인트
- Calculator는 사용자가 주입하는 구조 (ASE calculator interface)
  - 기본값으로 `xtb-python` 또는 외부 xTB ASE calculator 추천
  - 패키지 자체는 calculator-agnostic
- MD 엔진: ASE 내장 (`Langevin`, `VelocityVerlet`, `NPT`)
- 단계적 heating 옵션: `heating_schedule: list[tuple[float, int]]` (temp_K, steps)
- 로깅: trajectory 저장 옵션 (`equilibration_traj.traj`)
- PBC 설정: cell이 있으면 PBC 적용

---

## Step 4: Endpoint Swap & Relaxation

### 목적
평형화된 용매 환경을 유지하면서, 반응물/생성물 구조로 solute만 교체하고 interface를 relaxation한다.

### 상세

#### 4a. Geometry Swap
- 각 이미지 `reaction_images[i]`의 solute 좌표를 `equilibrated_system`의 solute 영역에 덮어씌움
- 용매(microsolv + bulk)는 그대로 유지

```python
# pseudo-code
def swap_solute(solvated_system, new_solute, n_solute):
    """Replace solute atoms in solvated system with new geometry."""
    new_system = solvated_system.copy()
    # Align new_solute to original solute position (center of mass alignment)
    # Replace positions[0:n_solute] with aligned new_solute positions
    new_system.positions[:n_solute] = align_to_reference(
        new_solute.positions,
        reference=solvated_system.positions[:n_solute]
    )
    return new_system
```

#### 4b. Constrained Relaxation
- Swap 후 solute-solvent interface에 clash가 있을 수 있으므로 relaxation 수행
- **방법 1**: Constrained optimization (solute 고정, 가까운 용매만 relax)
- **방법 2**: Short MD with soft constraints
- **방법 3**: FIRE/BFGS optimization with constraints

```python
# pseudo-code
solvated_images = []
for i, image in enumerate(reaction_images):
    swapped = swap_solute(equilibrated_system, image, n_solute)

    # Fix solute + distant solvent, relax interface
    # 또는: fix solute만 하고 전체 solvent relax
    swapped.set_constraint(FixAtoms(indices=list(range(n_solute))))
    swapped.calc = calculator

    if relax_method == "optimize":
        optimizer = BFGS(swapped)
        optimizer.run(fmax=fmax, steps=max_opt_steps)
    elif relax_method == "md_minimize":
        dyn = FIRE(swapped)
        dyn.run(fmax=fmax, steps=max_opt_steps)

    solvated_images.append(swapped.copy())
```

### 설계 포인트
- **Alignment**: Kabsch algorithm 또는 단순 center-of-mass alignment
  - 반응 경로에서 solute 구조가 크게 변하지 않는다고 가정
  - 큰 구조 변화 시 RMSD 기반 alignment (optional)
- **Relaxation strategy**: 3가지 중 선택 가능
  - `"optimize"`: BFGS/LBFGS (빠름, local minimum)
  - `"fire"`: FIRE optimizer (more robust)
  - `"md"`: short NVT MD (thermal relaxation)
- Solute는 항상 고정 (fix), 용매만 relax
- 중간 이미지들은 TS 기준 solvation이므로 swap + relax가 비교적 smooth

---

## Configuration Design

```python
@dataclass
class SolvationWorkflowConfig:
    """전체 워크플로우 설정"""

    # Step 1: Microsolvation
    microsolv_config: MicrosolvatorConfig  # 기존 microsolvator config 재사용
    ts_index: Optional[int] = None         # None이면 자동 (중간)

    # Step 2: Packmol
    box_margin: float = 10.0               # Å
    solvent_density: Optional[float] = None # g/cm³, None이면 n_solvent 사용
    n_bulk_solvent: Optional[int] = None    # 직접 지정 시
    min_distance: float = 2.0              # Å
    packmol_executable: Optional[str] = None

    # Step 3: MD Equilibration
    md_temperature: float = 300.0          # K
    md_steps: int = 5000
    md_timestep: float = 1.0              # fs
    md_friction: float = 0.01             # 1/fs (Langevin)
    heating_schedule: Optional[list] = None  # [(T, steps), ...]

    # Step 4: Relaxation
    relax_method: str = "optimize"         # "optimize" | "fire" | "md"
    relax_fmax: float = 0.05              # eV/Å
    relax_max_steps: int = 500

    # General
    working_directory: Optional[Path] = None
    keep_temps: bool = False
```

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
        calculator: Calculator,            # ASE calculator (pluggable)
        log_callback: Optional[Callable] = None,
    ) -> SolvatedTrajectoryResult:
        """Full pipeline: microsolvate → pack → equilibrate → swap & relax."""
        ...

    # 개별 step도 독립 호출 가능
    @staticmethod
    def microsolvate_ts(ts_image, solvent, config) -> MicrosolvationResult:
        """Step 1 only."""
        ...

    @staticmethod
    def pack_solvent_box(cluster, solvent, config) -> Atoms:
        """Step 2 only."""
        ...

    @staticmethod
    def equilibrate(system, n_fixed, calculator, config) -> Atoms:
        """Step 3 only."""
        ...

    @staticmethod
    def swap_and_relax(
        template, reaction_images, n_solute, calculator, config
    ) -> list[Atoms]:
        """Step 4 only."""
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
│   ├── config.py                 # SolvationWorkflowConfig
│   ├── builder.py                # SolvatedTrajectoryBuilder (main orchestrator)
│   ├── results.py                # SolvatedTrajectoryResult
│   ├── packmol.py                # Packmol wrapper (input gen, execution, parsing)
│   ├── equilibration.py          # MD equilibration logic
│   ├── swap.py                   # Geometry swap & alignment utilities
│   └── utils.py                  # Box size calculation, density estimation, etc.
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
| Step 3 (MD) | ASE + Calculator | calculator는 사용자 주입 |
| Step 4 (Swap/Relax) | ASE + Calculator | 동일 |
| Alignment | `scipy` (optional) | Kabsch rotation |

`pyproject.toml` 추가 의존성:
- `scipy` (optional, alignment용)
- calculator 관련은 사용자 책임 (e.g., `xtb-python`, `tblite` 등)

---

## Usage Example

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig
from microsolvator.workflow import (
    SolvatedTrajectoryBuilder,
    SolvationWorkflowConfig,
)

# 1. 반응 경로 로드
reaction_images = read("neb_initial_guess.traj", index=":")
solvent = read("water.xyz")

# 2. Calculator 준비 (사용자 선택)
from xtb.ase.calculator import XTB
calculator = XTB(method="GFN-FF")

# 3. 워크플로우 설정
config = SolvationWorkflowConfig(
    microsolv_config=MicrosolvatorConfig(
        nsolv=5,
        implicit_model="alpb",
        implicit_solvent="h2o",
        threads=8,
    ),
    box_margin=12.0,
    solvent_density=1.0,
    md_steps=10000,
    md_temperature=300.0,
    relax_method="optimize",
    relax_fmax=0.05,
)

# 4. 실행
result = SolvatedTrajectoryBuilder.build(
    reaction_images=reaction_images,
    solvent=solvent,
    config=config,
    calculator=calculator,
)

# 5. 결과
solvated_traj = result.solvated_images  # list[Atoms]
write("solvated_neb_guess.traj", solvated_traj)
```

---

## Edge Cases & 고려사항

1. **원자 순서 일관성**: 모든 step에서 `[solute | microsolv | bulk]` 순서 유지가 핵심
2. **Charged systems**: charge/multiplicity 전달 체계
3. **Large conformational change**: reactant↔product 구조 차이가 큰 경우 alignment 주의
4. **Packmol failure**: 밀도가 너무 높거나 box가 작으면 실패 → 자동 retry with larger box
5. **PBC consistency**: MD에서 PBC 사용 시 cell 정보가 모든 이미지에 전달되어야 함
6. **Calculator compatibility**: ASE calculator interface만 준수하면 어떤 calculator든 사용 가능
