"""Main orchestrator for the solvated reaction trajectory workflow."""

from __future__ import annotations

import copy
import tempfile
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, Union

from ase import Atoms

from ..results import MicrosolvationResult
from ..runner import Microsolvator
from .config import (
    EquilibrationConfig,
    MicrosolvatorConfig,
    PackmolConfig,
    RelaxationConfig,
    SolvationWorkflowConfig,
)
from .equilibration import equilibrate
from .packmol import PackmolSolvator
from .results import SolvatedTrajectoryResult
from .swap import relax_interface, swap_solute


# Type alias: either a calculator instance or a zero-arg callable that
# returns a fresh calculator.  Using a factory avoids sharing state
# across MD / optimization steps.
CalculatorLike = Any


def _make_calculator(calc: CalculatorLike) -> Any:
    """Return a fresh calculator instance.

    If *calc* is callable (factory), call it.  Otherwise deep-copy
    the instance so that each consumer gets independent state.
    """
    if callable(calc) and not hasattr(calc, "get_potential_energy"):
        # It's a factory function, not a calculator instance
        return calc()
    return copy.deepcopy(calc)


# ── Convenience function ────────────────────────────────────────────────────


def solvate_trajectory(
    reaction_images: List[Atoms],
    solvent: Atoms,
    *,
    calc: CalculatorLike,
    nsolv: int = 5,
    solvent_density: float = 1.0,
    ts_index: Optional[int] = None,
    config: Optional[SolvationWorkflowConfig] = None,
    **kwargs,
) -> SolvatedTrajectoryResult:
    """One-call interface for solvating a reaction trajectory.

    Parameters
    ----------
    reaction_images:
        Gas-phase NEB / string-method trajectory (list of ASE ``Atoms``).
    solvent:
        A single solvent molecule.
    calc:
        ASE calculator **instance** or a zero-argument **factory** that
        returns a fresh calculator.  Using a factory is recommended for
        stateful calculators (GPAW, VASP, …)::

            calc=lambda: XTB(method="GFN-FF")

    nsolv:
        Number of explicit microsolvation-shell solvent molecules
        (CREST ``--nsolv``).  Ignored when *config* is provided.
    solvent_density:
        Target bulk solvent density in g/cm³.  Ignored when *config*
        is provided.
    ts_index:
        TS image index.  ``None`` → ``len(images) // 2``.
    config:
        Full :class:`SolvationWorkflowConfig` for fine-grained control.
        When given, *nsolv* and *solvent_density* are ignored.
    **kwargs:
        Forwarded to :meth:`SolvatedTrajectoryBuilder.build`
        (e.g. ``log_callback``).

    Returns
    -------
    SolvatedTrajectoryResult

    Examples
    --------
    >>> from ase.io import read
    >>> from xtb.ase.calculator import XTB
    >>> from microsolvator.workflow import solvate_trajectory
    >>>
    >>> images = read("neb_guess.traj", index=":")
    >>> water  = read("water.xyz")
    >>> result = solvate_trajectory(
    ...     images, water,
    ...     calc=lambda: XTB(method="GFN-FF"),
    ...     nsolv=5,
    ... )
    """
    if config is None:
        config = SolvationWorkflowConfig(
            microsolv=MicrosolvatorConfig(nsolv=nsolv),
            packmol=PackmolConfig(solvent_density=solvent_density),
            ts_index=ts_index,
        )
    elif ts_index is not None:
        config.ts_index = ts_index

    return SolvatedTrajectoryBuilder.build(
        reaction_images=reaction_images,
        solvent=solvent,
        config=config,
        calc=calc,
        **kwargs,
    )


# ── Builder class ───────────────────────────────────────────────────────────


class SolvatedTrajectoryBuilder:
    """Build explicitly solvated reaction trajectories.

    The full pipeline is::

        microsolvate_ts  →  pack_solvent_box  →  equilibrate  →  swap_and_relax

    Each step is also exposed as an independent ``@staticmethod`` so that
    intermediate results can be reused or inspected without re-running the
    whole workflow.

    The *calc* parameter accepts either a calculator **instance** (which is
    deep-copied for each step) or a zero-argument **factory** callable.
    """

    @classmethod
    def build(
        cls,
        *,
        reaction_images: List[Atoms],
        solvent: Atoms,
        config: SolvationWorkflowConfig,
        calc: CalculatorLike = None,
        log_callback: Optional[Callable[[str, Any], None]] = None,
        # Backwards compatibility: accept 'calculator' as alias
        calculator: CalculatorLike = None,
    ) -> SolvatedTrajectoryResult:
        """Run the full solvation workflow.

        Parameters
        ----------
        reaction_images:
            Initial NEB / string-method trajectory.
        solvent:
            One solvent molecule.
        config:
            Full workflow configuration.
        calc:
            ASE calculator instance or factory callable.
        log_callback:
            Optional hook called at each step with
            ``(step_name, info_dict)``.
        """
        # Allow 'calculator' kwarg for backwards compat
        if calc is None and calculator is not None:
            calc = calculator
        if calc is None:
            raise TypeError(
                "A calculator is required. Pass calc=<calculator> or "
                "calc=lambda: <Calculator>(...)."
            )

        if not reaction_images:
            raise ValueError("reaction_images must contain at least one structure")

        ts_index = config.ts_index
        if ts_index is None:
            ts_index = len(reaction_images) // 2
        if not (0 <= ts_index < len(reaction_images)):
            raise ValueError(
                f"ts_index {ts_index} out of range "
                f"for {len(reaction_images)} images"
            )

        workdir = _resolve_workdir(config)

        # ── Step 1: Microsolvation ──────────────────────────────────────
        if log_callback:
            log_callback("microsolvation", {"ts_index": ts_index})

        ts_image = reaction_images[ts_index]
        microsolv_result = cls.microsolvate_ts(
            ts_image=ts_image,
            solvent=solvent,
            config=config.microsolv,
            working_directory=workdir / "microsolvation" if workdir else None,
        )

        microsolvated_ts = microsolv_result.best_structure or microsolv_result.final
        if microsolvated_ts is None:
            raise ValueError(
                "Microsolvation produced no valid cluster structure. "
                "Check CREST output in the working directory."
            )

        n_solute = len(ts_image)
        n_cluster = len(microsolvated_ts)

        # ── Step 2: Packmol box solvation ───────────────────────────────
        if log_callback:
            log_callback("packmol", {"n_cluster": n_cluster})

        boxed_system, _ = cls.pack_solvent_box(
            cluster=microsolvated_ts,
            solvent=solvent,
            config=config.packmol,
            workdir=workdir / "packmol" if workdir else None,
        )

        charge = ts_image.info.get("charge", config.microsolv.charge)
        multiplicity = ts_image.info.get("multiplicity", config.microsolv.uhf + 1)
        boxed_system.info["charge"] = charge
        boxed_system.info["multiplicity"] = multiplicity

        # ── Step 3: MD equilibration ────────────────────────────────────
        if log_callback:
            log_callback("equilibration", {"n_fixed": n_solute})

        traj_path = str(workdir / "equilibration.traj") if workdir else None
        equilibrated_system = cls.equilibrate(
            system=boxed_system,
            n_fixed=n_solute,
            calc=calc,
            config=config.equilibration,
            traj_path=traj_path,
        )

        # ── Step 4: Swap endpoints & relax ──────────────────────────────
        if log_callback:
            log_callback("swap_relax", {"n_images": len(reaction_images)})

        solvated_images = cls.swap_and_relax(
            template=equilibrated_system,
            reaction_images=reaction_images,
            n_solute=n_solute,
            ts_index=ts_index,
            calc=calc,
            config=config.relaxation,
        )

        for img in solvated_images:
            img.cell = equilibrated_system.cell.copy()
            img.pbc = equilibrated_system.pbc.copy()
            img.info["charge"] = charge
            img.info["multiplicity"] = multiplicity

        return SolvatedTrajectoryResult(
            solvated_images=solvated_images,
            microsolvated_ts=microsolvated_ts,
            boxed_system=boxed_system,
            equilibrated_system=equilibrated_system,
            n_solute_atoms=n_solute,
            n_cluster_atoms=n_cluster,
            n_total_atoms=len(equilibrated_system),
            ts_index=ts_index,
            config=config,
            working_directory=workdir or Path("."),
        )

    # ── Individual step methods ─────────────────────────────────────────

    @staticmethod
    def microsolvate_ts(
        ts_image: Atoms,
        solvent: Atoms,
        config: MicrosolvatorConfig,
        **kwargs,
    ) -> MicrosolvationResult:
        """Step 1: CREST QCG microsolvation of the TS image."""
        return Microsolvator.run(
            solute=ts_image,
            solvent=solvent,
            config=config,
            constrain_solute=True,
            **kwargs,
        )

    @staticmethod
    def pack_solvent_box(
        cluster: Atoms,
        solvent: Atoms,
        config: PackmolConfig,
        workdir: Optional[Path] = None,
    ) -> Tuple[Atoms, int]:
        """Step 2: Packmol box solvation."""
        return PackmolSolvator.run(
            cluster=cluster,
            solvent=solvent,
            config=config,
            workdir=workdir,
        )

    @staticmethod
    def equilibrate(
        system: Atoms,
        n_fixed: int,
        calc: CalculatorLike = None,
        config: Optional[EquilibrationConfig] = None,
        traj_path: Optional[str] = None,
        log_callback: Optional[Callable[[str, Any], None]] = None,
        *,
        calculator: CalculatorLike = None,
    ) -> Atoms:
        """Step 3: Constrained MD equilibration.

        A fresh calculator is created from *calc* to avoid state leakage.
        """
        c = calc if calc is not None else calculator
        if config is None:
            config = EquilibrationConfig()
        return equilibrate(
            system=system,
            n_fixed=n_fixed,
            calculator=_make_calculator(c),
            config=config,
            traj_path=traj_path,
            log_callback=log_callback,
        )

    @staticmethod
    def swap_and_relax(
        template: Atoms,
        reaction_images: List[Atoms],
        n_solute: int,
        ts_index: int,
        calc: CalculatorLike = None,
        config: Optional[RelaxationConfig] = None,
        *,
        calculator: CalculatorLike = None,
    ) -> List[Atoms]:
        """Step 4: Swap each image's solute and relax the solvent.

        Propagates outward from the TS image in both directions.  Each
        image gets a **fresh calculator copy** to avoid state leakage.
        """
        c = calc if calc is not None else calculator
        if config is None:
            config = RelaxationConfig()

        n_images = len(reaction_images)
        results: List[Optional[Atoms]] = [None] * n_images

        results[ts_index] = template.copy()

        # Forward: TS+1, TS+2, … → product
        current = template
        for i in range(ts_index + 1, n_images):
            swapped = swap_solute(current, reaction_images[i], n_solute)
            relaxed = relax_interface(
                swapped, n_solute, _make_calculator(c), config
            )
            results[i] = relaxed
            current = relaxed

        # Backward: TS-1, TS-2, … → reactant
        current = template
        for i in range(ts_index - 1, -1, -1):
            swapped = swap_solute(current, reaction_images[i], n_solute)
            relaxed = relax_interface(
                swapped, n_solute, _make_calculator(c), config
            )
            results[i] = relaxed
            current = relaxed

        return results  # type: ignore[return-value]


def _resolve_workdir(config: SolvationWorkflowConfig) -> Optional[Path]:
    """Return a working directory Path, or None for fully in-memory execution."""
    if config.working_directory is not None:
        p = Path(config.working_directory)
        p.mkdir(parents=True, exist_ok=True)
        return p
    if config.keep_temps:
        p = Path(tempfile.mkdtemp(prefix="solvated_traj_"))
        return p
    return None
