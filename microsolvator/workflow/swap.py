"""Geometry swap and interface relaxation utilities (Step 4)."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin
import ase.units as units

from .config import RelaxationConfig


def kabsch_align(mobile: Atoms, target: Atoms) -> Atoms:
    """Return a copy of *mobile* Kabsch-aligned onto *target*.

    Minimises the RMSD between *mobile* and *target* by finding the optimal
    rotation matrix (Kabsch algorithm) followed by a translation to match
    the centres of mass.

    Parameters
    ----------
    mobile:
        Structure to be rotated / translated.
    target:
        Reference structure (not modified).

    Returns
    -------
    Atoms
        Copy of *mobile* with positions replaced by the aligned coordinates.

    Raises
    ------
    ValueError
        If atom counts or element lists differ between *mobile* and *target*.
    """
    if len(mobile) != len(target):
        raise ValueError(
            f"Cannot align structures with different atom counts: "
            f"{len(mobile)} vs {len(target)}"
        )
    if list(mobile.symbols) != list(target.symbols):
        raise ValueError(
            "Atom species mismatch between mobile and target.\n"
            f"mobile : {list(mobile.symbols)[:8]}…\n"
            f"target : {list(target.symbols)[:8]}…"
        )

    P = mobile.positions.copy()   # (N, 3) – mobile
    Q = target.positions.copy()   # (N, 3) – target

    p_com = P.mean(axis=0)
    q_com = Q.mean(axis=0)

    P_c = P - p_com
    Q_c = Q - q_com

    # Covariance matrix
    H = P_c.T @ Q_c

    # SVD
    U, _, Vt = np.linalg.svd(H)

    # Correct for improper rotation (reflection)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T

    aligned = mobile.copy()
    aligned.positions = (P_c @ R) + q_com
    return aligned


def swap_solute(
    solvated_system: Atoms,
    new_solute: Atoms,
    n_solute: int,
) -> Atoms:
    """Replace the solute in *solvated_system* with a Kabsch-aligned *new_solute*.

    The solvent atoms (indices ``n_solute`` onward) are left untouched.

    Parameters
    ----------
    solvated_system:
        Equilibrated solvated box used as the template.
    new_solute:
        New solute geometry (e.g., a reactant or product image).
    n_solute:
        Number of solute atoms at the front of *solvated_system*.

    Returns
    -------
    Atoms
        New system with *new_solute* aligned into the solute slot.
    """
    ref_solute = solvated_system[:n_solute]
    aligned = kabsch_align(new_solute, ref_solute)

    new_system = solvated_system.copy()
    new_system.positions[:n_solute] = aligned.positions
    return new_system


def relax_interface(
    system: Atoms,
    n_solute: int,
    calculator,
    config: RelaxationConfig,
) -> Atoms:
    """Fix the solute and relax the solvent to remove steric clashes.

    Parameters
    ----------
    system:
        Swapped system (solute positions updated, solvent positions from template).
    n_solute:
        Number of leading atoms to freeze during relaxation.
    calculator:
        ASE-compatible calculator.
    config:
        Relaxation settings (method, convergence criteria, etc.).

    Returns
    -------
    Atoms
        Relaxed copy of *system*.
    """
    system = system.copy()
    system.set_constraint(FixAtoms(indices=list(range(n_solute))))
    system.calc = calculator

    method = config.method.lower()

    if method in ("optimize", "lbfgs"):
        from ase.optimize import LBFGS

        opt = LBFGS(system, logfile=None)
        opt.run(fmax=config.fmax, steps=config.max_steps)

    elif method == "fire":
        from ase.optimize import FIRE

        opt = FIRE(system, logfile=None)
        opt.run(fmax=config.fmax, steps=config.max_steps)

    elif method == "md":
        dyn = Langevin(
            system,
            timestep=config.md_timestep * units.fs,
            temperature_K=config.md_temperature,
            friction=config.md_friction,
        )
        dyn.run(config.md_steps)

    else:
        raise ValueError(
            f"Unknown relaxation method {config.method!r}. "
            "Choose 'optimize', 'fire', or 'md'."
        )

    return system
