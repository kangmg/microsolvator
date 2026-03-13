"""Configuration dataclasses for the solvation workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from ..config import MicrosolvatorConfig


@dataclass
class PackmolConfig:
    """Configuration for Packmol box solvation (Step 2).

    Provide exactly one of *solvent_density* or *n_bulk_solvent*.
    Validation is deferred to runtime so that the dataclass can be
    constructed without either value and filled in later.
    """

    box_margin: float = 10.0
    """Å of padding from the cluster bounding box to the box wall."""

    solvent_density: Optional[float] = None
    """Target bulk solvent density in g/cm³ (mutually exclusive with n_bulk_solvent)."""

    n_bulk_solvent: Optional[int] = None
    """Explicit number of bulk solvent molecules (mutually exclusive with solvent_density)."""

    min_distance: float = 2.0
    """Minimum distance (Å) between any two molecules in the Packmol output."""

    packmol_executable: Optional[str] = None
    """Path to the packmol binary; None means search PATH."""

    retry_margin_increment: float = 2.0
    """Extra Å added to box_margin on each Packmol retry."""

    max_retries: int = 3
    """Maximum number of Packmol retries after failure."""


@dataclass
class EquilibrationConfig:
    """Configuration for constrained MD equilibration (Step 3).

    Default schedule: NVT heating (if heating_schedule provided) →
    NPT density relaxation (if npt_steps > 0) → NVT production.
    """

    temperature: float = 300.0
    """Target temperature in K for NVT production (and NPT)."""

    nvt_steps: int = 5000
    """Number of NVT production steps."""

    npt_steps: int = 0
    """Number of NPT density-relaxation steps. 0 disables the NPT phase."""

    npt_pfactor: Optional[float] = None
    """NPT barostat coupling factor (eV·fs²/Å³).
    None → auto-estimated from a 2.2 GPa bulk modulus (water-like)."""

    timestep: float = 1.0
    """MD timestep in fs."""

    friction: float = 0.01
    """Langevin friction coefficient in 1/fs."""

    heating_schedule: Optional[List[Tuple[float, int]]] = None
    """Stepwise NVT heating: list of (temperature_K, n_steps) pairs.
    Executed in order before the NPT / NVT production phases."""


@dataclass
class RelaxationConfig:
    """Configuration for endpoint swap & interface relaxation (Step 4)."""

    method: str = "optimize"
    """Relaxation method: 'optimize' (LBFGS), 'fire', or 'md' (Langevin)."""

    fmax: float = 0.05
    """Force convergence threshold in eV/Å (used by 'optimize' and 'fire')."""

    max_steps: int = 500
    """Maximum optimizer steps (used by 'optimize' and 'fire')."""

    md_steps: int = 200
    """Number of MD steps (used when method='md')."""

    md_temperature: float = 300.0
    """MD temperature in K (used when method='md')."""

    md_timestep: float = 0.5
    """MD timestep in fs (used when method='md')."""

    md_friction: float = 0.01
    """Langevin friction in 1/fs (used when method='md')."""


@dataclass
class SolvationWorkflowConfig:
    """Top-level configuration that composes all step-specific configs.

    Example
    -------
    >>> config = SolvationWorkflowConfig(
    ...     microsolv=MicrosolvatorConfig(nsolv=5),
    ...     packmol=PackmolConfig(solvent_density=1.0),
    ... )
    """

    microsolv: MicrosolvatorConfig
    """CREST microsolvation config (Step 1). Controls ensemble flag, nsolv, etc."""

    packmol: PackmolConfig = field(default_factory=PackmolConfig)
    """Packmol box solvation config (Step 2)."""

    equilibration: EquilibrationConfig = field(default_factory=EquilibrationConfig)
    """MD equilibration config (Step 3)."""

    relaxation: RelaxationConfig = field(default_factory=RelaxationConfig)
    """Endpoint relaxation config (Step 4)."""

    ts_index: Optional[int] = None
    """Index of the TS image in reaction_images. None → len(images) // 2."""

    working_directory: Optional[Path] = None
    """Directory for intermediate files. None → system temp (removed after run)."""

    keep_temps: bool = False
    """Keep intermediate files when working_directory is None."""
