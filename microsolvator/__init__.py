"""Public API for the microsolvator package."""

__version__ = "0.1.0"

from .config import MicrosolvatorConfig
from .results import MicrosolvationResult
from .runner import Microsolvator
from .support import (
    list_supported_implicit_solvents,
    supports_implicit_solvent,
)
from .install import (
    PACKAGE_BIN_DIR,
    install_crest,
    install_xtb,
    resolve_crest_binary,
    resolve_xtb_binary,
)

__all__ = [
    "__version__",
    "MicrosolvatorConfig",
    "MicrosolvationResult",
    "Microsolvator",
    "list_supported_implicit_solvents",
    "supports_implicit_solvent",
    "install_crest",
    "install_xtb",
    "resolve_crest_binary",
    "resolve_xtb_binary",
    "PACKAGE_BIN_DIR",
]
