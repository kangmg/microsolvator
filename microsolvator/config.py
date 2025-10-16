"""Configuration objects for microsolvator runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence


@dataclass
class MicrosolvatorConfig:
    """User-facing configuration for a CREST microsolvation run."""

    nsolv: int
    method: str = "gfn2"
    temperature: float = 298.0
    mdtime: float = 50.0
    charge: int = 0
    uhf: int = 0
    crest_executable: Optional[str] = None
    xtb_executable: Optional[str] = None
    implicit_model: Optional[str] = None
    implicit_solvent: Optional[str] = None
    ensemble: bool = True
    nopreopt: bool = False
    additional_flags: Sequence[str] = field(default_factory=tuple)

    def build_flag_list(self) -> List[str]:
        """Return a mutable list of CREST CLI flags without solute/solvent files."""

        flags: List[str] = ["--qcg"]
        if self.nsolv <= 0:
            raise ValueError("nsolv must be a positive integer")
        flags.extend(["--nsolv", str(self.nsolv)])

        flags.extend(["--T", "%g" % self.temperature])
        flags.extend(["--mdtime", str(self.mdtime)])

        flags.extend(_method_flags(self.method))

        if self.ensemble:
            flags.append("--ensemble")

        if self.implicit_model and self.implicit_solvent:
            flags.extend([f"--{self.implicit_model}", self.implicit_solvent])
        elif self.implicit_model or self.implicit_solvent:
            raise ValueError("Both implicit_model and implicit_solvent must be provided together")

        flags.extend(["--chrg", str(self.charge)])
        flags.extend(["--uhf", str(self.uhf)])

        flags.extend(self.additional_flags)

        if self.nopreopt:
            flags.append("--nopreopt")

        return flags

    @classmethod
    def from_kwargs(cls, **kwargs: object) -> "MicrosolvatorConfig":
        """Helper for creating configs from keyword arguments."""

        return cls(**kwargs)  # type: ignore[arg-type]


def _method_flags(method: str) -> List[str]:
    name = method.strip().lower()

    if name in {"gfn2", "gfn-2"}:
        return ["--gfn", "2"]
    if name in {"gfn1", "gfn-1"}:
        return ["--gfn", "1"]
    if name in {"gfn0", "gfn-0"}:
        return ["--gfn", "0"]
    if name in {"gfnff", "gfn-ff", "ff"}:
        return ["--gfnff"]

    raise ValueError(f"Unsupported method '{method}' for CREST run")
