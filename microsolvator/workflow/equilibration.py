"""Constrained MD equilibration (Step 3)."""

from __future__ import annotations

from typing import Callable, Optional

import ase.units as units
from ase import Atoms
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin

from .config import EquilibrationConfig


def equilibrate(
    system: Atoms,
    n_fixed: int,
    calculator,
    config: EquilibrationConfig,
    traj_path: Optional[str] = None,
    log_callback: Optional[Callable[[str, object], None]] = None,
) -> Atoms:
    """Run constrained MD to equilibrate the bulk solvent.

    Atoms ``0 … n_fixed-1`` (the microsolvated cluster) are frozen via
    :class:`~ase.constraints.FixAtoms`.  The remaining bulk solvent atoms
    are free to move.

    Equilibration proceeds in up to three phases:

    1. **NVT heating** — optional stepwise ramp defined by
       ``config.heating_schedule``.
    2. **NPT density relaxation** — optional; runs when
       ``config.npt_steps > 0``.
    3. **NVT production** — always runs for ``config.nvt_steps`` steps.

    Parameters
    ----------
    system:
        The Packmol-generated solvated box (modified in-place on a copy).
    n_fixed:
        Number of leading atoms to freeze (solute + microsolvation shell).
    calculator:
        Any ASE-compatible calculator.  Charge / multiplicity should be
        configured on the calculator before passing.
    config:
        Equilibration settings.
    traj_path:
        If given, write an ASE Trajectory to this path (every 10 steps).
    log_callback:
        Optional hook called every 10 steps with ``("md", atoms)``.

    Returns
    -------
    Atoms
        A copy of the equilibrated system (constraints removed from copy).
    """
    system = system.copy()
    system.set_constraint(FixAtoms(indices=list(range(n_fixed))))
    system.calc = calculator

    if system.cell.any() and not system.pbc.any():
        system.pbc = True

    traj_writer = None
    if traj_path:
        from ase.io.trajectory import Trajectory

        traj_writer = Trajectory(traj_path, "w", system)

    def _attach(dyn: object) -> None:
        if traj_writer is not None:
            dyn.attach(traj_writer.write, interval=10)  # type: ignore[attr-defined]
        if log_callback is not None:
            dyn.attach(  # type: ignore[attr-defined]
                lambda: log_callback("md", system), interval=10
            )

    try:
        # Phase 1: NVT heating
        if config.heating_schedule:
            for temp_K, steps in config.heating_schedule:
                dyn = Langevin(
                    system,
                    timestep=config.timestep * units.fs,
                    temperature_K=temp_K,
                    friction=config.friction,
                )
                _attach(dyn)
                dyn.run(steps)

        # Phase 2: NPT density relaxation (optional)
        if config.npt_steps > 0:
            from ase.md.npt import NPT

            pfactor = config.npt_pfactor
            if pfactor is None:
                # Default: ~2.2 GPa bulk modulus (water-like), ttime = 25 fs
                ttime = 25.0 * units.fs
                pfactor = 2.2 * units.GPa * ttime ** 2

            dyn = NPT(
                system,
                timestep=config.timestep * units.fs,
                temperature_K=config.temperature,
                externalstress=0.0,
                ttime=25.0 * units.fs,
                pfactor=pfactor,
            )
            _attach(dyn)
            dyn.run(config.npt_steps)

        # Phase 3: NVT production
        dyn = Langevin(
            system,
            timestep=config.timestep * units.fs,
            temperature_K=config.temperature,
            friction=config.friction,
        )
        _attach(dyn)
        dyn.run(config.nvt_steps)

    finally:
        if traj_writer is not None:
            traj_writer.close()

    return system
