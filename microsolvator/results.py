"""Result containers for microsolvation runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

from ase import Atoms


@dataclass
class MicrosolvationResult:
    """Structured output for a CREST microsolvation calculation."""

    command: Sequence[str]
    working_directory: Path
    best_structure: Optional[Atoms]
    ensemble: List[Atoms]
    population_path: Optional[Path]
    stdout: str
    stderr: str

    def ensure_outputs(self) -> None:
        """Raise if no ensemble was parsed."""

        if not self.ensemble and self.best_structure is None:
            raise ValueError("Microsolvation run produced no structures")
