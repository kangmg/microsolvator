"""Result container for the solvation workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

from ase import Atoms

from .config import SolvationWorkflowConfig


@dataclass
class SolvatedTrajectoryResult:
    """Complete output of the solvation workflow.

    Attributes
    ----------
    solvated_images:
        Final explicitly solvated trajectory.  Atom ordering in every image
        is ``[solute | microsolvation_shell | bulk_solvent]``.
    microsolvated_ts:
        CREST QCG cluster from Step 1 (TS + first solvation shell).
    boxed_system:
        Packmol output from Step 2 (cluster + bulk solvent, with cell).
    equilibrated_system:
        MD-equilibrated system from Step 3 (template for Step 4 swaps).
    n_solute_atoms:
        Number of atoms belonging to the solute (= len of one NEB image).
    n_cluster_atoms:
        ``n_solute_atoms`` + microsolvation shell atoms.
    n_total_atoms:
        Total atoms in the solvated box.
    ts_index:
        The TS image index that was actually used (resolved from config).
    config:
        The workflow configuration used for this run.
    working_directory:
        Directory where intermediate files were written.
    """

    solvated_images: List[Atoms]
    microsolvated_ts: Atoms
    boxed_system: Atoms
    equilibrated_system: Atoms
    n_solute_atoms: int
    n_cluster_atoms: int
    n_total_atoms: int
    ts_index: int
    config: SolvationWorkflowConfig
    working_directory: Path
