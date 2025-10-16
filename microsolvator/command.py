"""Utilities to construct CREST command invocations."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

from .config import MicrosolvatorConfig


def build_crest_command(
    *,
    config: MicrosolvatorConfig,
    solute_path: Path,
    solvent_path: Path,
) -> List[str]:
    """Return the full CREST command as a list suitable for subprocess calls."""

    command: List[str] = [config.crest_executable, str(solute_path)]
    for flag in config.build_flag_list():
        command.append(flag)
        if flag == "--qcg":
            command.append(str(solvent_path))
    return command
