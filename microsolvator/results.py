"""Result containers for microsolvation runs."""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field
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
    executed: bool = True
    final: Optional[Atoms] = None
    traj: List[Atoms] = field(default_factory=list)

    def ensure_outputs(self) -> None:
        """Raise if no ensemble was parsed for an executed run."""

        if not self.executed:
            return

        if not self.ensemble and self.best_structure is None and self.final is None:
            raise ValueError("Microsolvation run produced no structures")

    @property
    def shell_command(self) -> str:
        """Return the command as a shell-escaped string."""

        return shlex.join(self.command)
