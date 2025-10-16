"""Binary resolution and installation helpers for CREST and xTB."""

from __future__ import annotations

import os
import shutil
import stat
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional
from shutil import which


DEFAULT_CREST_URL = (
    "https://github.com/crest-lab/crest/releases/download/latest/"
    "crest-gnu-12-ubuntu-latest.tar.xz"
)
DEFAULT_XTB_URL = (
    "https://github.com/grimme-lab/xtb/releases/download/v6.7.0/"
    "xtb-6.7.0-linux-x86_64.tar.xz"
)

PACKAGE_BIN_DIR = Path(__file__).resolve().parent / "_bin"

__all__ = [
    "PACKAGE_BIN_DIR",
    "install_crest",
    "install_xtb",
    "resolve_crest_binary",
    "resolve_xtb_binary",
]


def resolve_crest_binary(explicit_path: Optional[str]) -> str:
    """Resolve the CREST executable path following the configured priority order."""

    return _resolve_binary(
        explicit_path=explicit_path,
        env_var="CREST_BIN",
        binary_name="crest",
    )


def resolve_xtb_binary(explicit_path: Optional[str]) -> str:
    """Resolve the xTB executable path following the configured priority order."""

    return _resolve_binary(
        explicit_path=explicit_path,
        env_var="XTB_BIN",
        binary_name="xtb",
    )


def install_crest(
    *, url: str = DEFAULT_CREST_URL, force: bool = False
) -> Path:
    """Download and install CREST into the package bin directory."""

    return _install_binary(
        url=url,
        archive_name="crest.tar.xz",
        expected_binary="crest",
        force=force,
    )


def install_xtb(
    *, url: str = DEFAULT_XTB_URL, force: bool = False
) -> Path:
    """Download and install xTB into the package bin directory."""

    return _install_binary(
        url=url,
        archive_name="xtb.tar.xz",
        expected_binary="xtb",
        force=force,
    )


def _resolve_binary(*, explicit_path: Optional[str], env_var: str, binary_name: str) -> str:
    if explicit_path:
        explicit = Path(explicit_path).expanduser()
        if not explicit.exists():
            raise FileNotFoundError(f"Specified {binary_name} executable not found: {explicit}")
        return str(explicit.resolve())

    env_path = os.environ.get(env_var)
    if env_path:
        env_candidate = Path(env_path).expanduser()
        if not env_candidate.exists():
            raise FileNotFoundError(
                f"Environment variable {env_var} points to missing {binary_name}: {env_candidate}"
            )
        return str(env_candidate.resolve())

    package_candidate = _find_in_package(binary_name)
    if package_candidate is not None:
        return str(package_candidate.resolve())

    located = which(binary_name)
    if located:
        return str(Path(located).resolve())

    return binary_name


def _find_in_package(binary_name: str) -> Optional[Path]:
    if not PACKAGE_BIN_DIR.exists():
        return None

    for candidate in PACKAGE_BIN_DIR.rglob(binary_name):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return None


def _install_binary(
    *,
    url: str,
    archive_name: str,
    expected_binary: str,
    force: bool,
) -> Path:
    PACKAGE_BIN_DIR.mkdir(parents=True, exist_ok=True)
    archive_path = PACKAGE_BIN_DIR / archive_name

    if archive_path.exists():
        archive_path.unlink()

    urllib.request.urlretrieve(url, archive_path)

    extract_dir = PACKAGE_BIN_DIR / f"_{expected_binary}_extract"
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir()

    with tarfile.open(archive_path, mode="r:xz") as tar:
        tar.extractall(path=extract_dir)

    archive_path.unlink(missing_ok=True)

    binary_path = _promote_binary(extract_dir, expected_binary)

    if force and (PACKAGE_BIN_DIR / expected_binary).exists():
        (PACKAGE_BIN_DIR / expected_binary).unlink()

    if binary_path is None:
        raise FileNotFoundError(f"Failed to locate {expected_binary} inside extracted archive")

    final_path = PACKAGE_BIN_DIR / expected_binary
    if final_path.exists():
        final_path.unlink()
    shutil.move(str(binary_path), final_path)

    _ensure_executable(final_path)

    shutil.rmtree(extract_dir)

    return final_path


def _promote_binary(root: Path, binary_name: str) -> Optional[Path]:
    for candidate in root.rglob(binary_name):
        if candidate.is_file():
            return candidate
    return None


def _ensure_executable(path: Path) -> None:
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
