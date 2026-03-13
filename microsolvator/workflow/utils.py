"""Utility helpers for the solvation workflow."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import atomic_masses


# Avogadro's number (mol⁻¹)
_N_A = 6.02214076e23

# Å³ → cm³ conversion
_ANG3_TO_CM3 = 1e-24


def estimate_box_size(cluster: Atoms, margin: float) -> float:
    """Return the cubic box side length (Å) that fits *cluster* with *margin*.

    The span of the cluster along each Cartesian axis is computed; the
    largest span plus twice the margin gives the box side length.
    """
    pos = cluster.positions
    span = pos.max(axis=0) - pos.min(axis=0)
    return float(span.max()) + 2.0 * margin


def count_solvent_molecules(
    box_size: float,
    solvent: Atoms,
    density: float,
) -> int:
    """Estimate the number of solvent molecules to fill a cubic box.

    Parameters
    ----------
    box_size:
        Cubic box side length in Å.
    solvent:
        One solvent molecule (used to obtain its molar mass).
    density:
        Target density in g/cm³.

    Returns
    -------
    int
        Number of solvent molecules (at least 1).
    """
    vol_cm3 = (box_size ** 3) * _ANG3_TO_CM3
    molar_mass = float(sum(atomic_masses[n] for n in solvent.numbers))
    n = int(vol_cm3 * density * _N_A / molar_mass)
    return max(n, 1)


def validate_packmol_output(
    output: Atoms,
    cluster: Atoms,
    n_bulk_solvent: int,
    n_solvent_atoms: int,
) -> None:
    """Assert that Packmol output preserves the cluster at index 0.

    Checks both the total atom count and that the leading symbols match
    the cluster exactly.

    Raises
    ------
    ValueError
        If atom count or element ordering does not match expectations.
    """
    n_cluster = len(cluster)
    expected = n_cluster + n_bulk_solvent * n_solvent_atoms

    if len(output) != expected:
        raise ValueError(
            f"Packmol output atom count mismatch: "
            f"expected {expected} ({n_cluster} cluster + "
            f"{n_bulk_solvent}×{n_solvent_atoms} solvent), "
            f"got {len(output)}"
        )

    cluster_syms = list(cluster.symbols)
    output_syms = list(output.symbols[:n_cluster])
    if output_syms != cluster_syms:
        raise ValueError(
            "Atom ordering in Packmol output does not match cluster.\n"
            f"Expected first {n_cluster} symbols: {cluster_syms[:8]}…\n"
            f"Got:                                {output_syms[:8]}…"
        )
