"""Solvated reaction trajectory workflow.

Build explicitly solvated NEB / string-method trajectory guesses from an
initial gas-phase reaction path.

Typical usage
-------------
>>> from microsolvator.workflow import solvate_trajectory
>>> result = solvate_trajectory(
...     images, water,
...     calc=lambda: XTB(method="GFN-FF"),
...     nsolv=5,
... )
"""

from .builder import SolvatedTrajectoryBuilder, solvate_trajectory
from .config import (
    EquilibrationConfig,
    PackmolConfig,
    RelaxationConfig,
    SolvationWorkflowConfig,
)
from .results import SolvatedTrajectoryResult

__all__ = [
    "solvate_trajectory",
    "SolvatedTrajectoryBuilder",
    "SolvationWorkflowConfig",
    "PackmolConfig",
    "EquilibrationConfig",
    "RelaxationConfig",
    "SolvatedTrajectoryResult",
]
