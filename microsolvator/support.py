"""Implicit solvent support tables and helpers."""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple


_SOLVENT_TABLE: Dict[str, Dict[str, Dict[str, bool]]] = {
    "gfn1": {
        "alpb": {
            "acetone": True,
            "acetonitrile": True,
            "aniline": True,
            "benzaldehyde": True,
            "benzene": True,
            "ch2cl2": True,
            "chcl3": True,
            "cs2": True,
            "dioxane": True,
            "dmf": True,
            "dmso": True,
            "ether": True,
            "ethylacetate": True,
            "furane": True,
            "hexadecane": True,
            "hexane": True,
            "methanol": True,
            "nitromethane": True,
            "octanol": True,
            "phenol": True,
            "toluene": True,
            "thf": True,
            "h2o": True,
        },
        "gbsa": {
            "acetone": True,
            "acetonitrile": True,
            "aniline": False,
            "benzaldehyde": False,
            "benzene": True,
            "ch2cl2": True,
            "chcl3": True,
            "cs2": True,
            "dioxane": False,
            "dmf": False,
            "dmso": True,
            "ether": True,
            "ethylacetate": False,
            "furane": False,
            "hexadecane": False,
            "hexane": False,
            "methanol": True,
            "nitromethane": False,
            "octanol": False,
            "phenol": False,
            "toluene": True,
            "thf": True,
            "h2o": True,
        },
    },
    "gfn2": {
        "alpb": {
            "acetone": True,
            "acetonitrile": True,
            "aniline": True,
            "benzaldehyde": True,
            "benzene": True,
            "ch2cl2": True,
            "chcl3": True,
            "cs2": True,
            "dioxane": True,
            "dmf": True,
            "dmso": True,
            "ether": True,
            "ethylacetate": True,
            "furane": True,
            "hexadecane": True,
            "hexane": True,
            "methanol": True,
            "nitromethane": True,
            "octanol": True,
            "phenol": True,
            "toluene": True,
            "thf": True,
            "h2o": True,
        },
        "gbsa": {
            "acetone": True,
            "acetonitrile": True,
            "aniline": False,
            "benzaldehyde": False,
            "benzene": True,
            "ch2cl2": True,
            "chcl3": True,
            "cs2": True,
            "dioxane": False,
            "dmf": True,
            "dmso": True,
            "ether": True,
            "ethylacetate": False,
            "furane": False,
            "hexadecane": False,
            "hexane": True,
            "methanol": True,
            "nitromethane": False,
            "octanol": False,
            "phenol": False,
            "toluene": True,
            "thf": True,
            "h2o": True,
        },
    },
    "gfn-ff": {
        "alpb": {
            "acetone": True,
            "acetonitrile": True,
            "aniline": True,
            "benzaldehyde": True,
            "benzene": True,
            "ch2cl2": True,
            "chcl3": True,
            "cs2": True,
            "dioxane": True,
            "dmf": True,
            "dmso": True,
            "ether": True,
            "ethylacetate": True,
            "furane": True,
            "hexadecane": True,
            "hexane": True,
            "methanol": True,
            "nitromethane": True,
            "octanol": True,
            "phenol": True,
            "toluene": True,
            "thf": True,
            "h2o": True,
        }
    },
}


def supports_implicit_solvent(*, method: str, model: str, solvent: str) -> bool:
    """Return True if the method/model/solvent combination is supported."""

    method_key = method.lower()
    model_key = model.lower()
    solvent_key = solvent.lower()

    method_entry = _SOLVENT_TABLE.get(method_key)
    if method_entry is None:
        return False

    model_entry = method_entry.get(model_key)
    if model_entry is None:
        return False

    return model_entry.get(solvent_key, False)


def list_supported_implicit_solvents(
    *, method: Optional[str] = None, model: Optional[str] = None
) -> Dict[str, Dict[str, Tuple[str, ...]]]:
    """Return a filtered view of supported solvents.

    The result groups solvents by method and model where the support flag is True.
    """

    data: Dict[str, Dict[str, Tuple[str, ...]]] = {}
    for method_key, models in _SOLVENT_TABLE.items():
        if method and method_key != method.lower():
            continue
        data[method_key] = {}
        for model_key, solvents in models.items():
            if model and model_key != model.lower():
                continue
            supported = tuple(sorted(solvent for solvent, ok in solvents.items() if ok))
            if supported:
                data[method_key][model_key] = supported
        if not data[method_key]:
            data.pop(method_key)
    return data


def validate_implicit_choice(*, method: str, model: str, solvent: str) -> None:
    """Raise ValueError if the combination is unsupported."""

    if not supports_implicit_solvent(method=method, model=model, solvent=solvent):
        raise ValueError(
            "Unsupported implicit solvent combination: method=%s, model=%s, solvent=%s"
            % (method, model, solvent)
        )
