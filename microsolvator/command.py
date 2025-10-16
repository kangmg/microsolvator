"""Utilities to construct CREST command invocations."""

from __future__ import annotations

import os
from pathlib import Path
from typing import List
from shutil import which

from .config import MicrosolvatorConfig


def build_crest_command(
    *,
    config: MicrosolvatorConfig,
    crest_executable: str,
    xtb_executable: str,
    solute_path: Path,
    solvent_path: Path,
) -> List[str]:
    """Return the full CREST command as a list suitable for subprocess calls."""

    command: List[str] = [
        _normalize_exec_path(crest_executable),
        str(solute_path.resolve()),
    ]
    for flag in config.build_flag_list():
        command.append(flag)
        if flag == "--qcg":
            command.append(str(solvent_path.resolve()))
    command.extend(["--xnam", _normalize_exec_path(xtb_executable)])
    return command


def _normalize_exec_path(executable: str) -> str:
    path = Path(executable).expanduser()
    if path.is_absolute():
        return str(path.resolve())
    string_path = str(path)
    if os.sep in string_path or (os.altsep and os.altsep in string_path):
        return str(path.resolve())
    resolved = which(executable)
    if resolved:
        return str(Path(resolved).resolve())
    return executable
