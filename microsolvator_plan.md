% microsolvator Implementation Notes

## Implementation Plan

1. **Package Skeleton**
   - Create a Python package `microsolvator` with submodules for command construction, execution, and result parsing.
   - Expose a high-level API (e.g., `Microsolvator.run(...)`) that accepts ASE `Atoms` objects for solute and solvent plus run parameters (`nsolv`, `temperature`, `mdtime`, etc.).
   - Provide lightweight data classes for configuration (paths, crest executable options) and constraint specification.

2. **CREST Command Builder**
   - Serialize provided ASE `Atoms` to temporary XYZ files and generate `.xcontrol` files on-the-fly from supplied atom index sets or pre-built templates.
   - Provide helpers for global solute constraints by writing all solute atom indices into the generated `.xcontrol` file when requested.
   - Construct CREST command line arguments from configuration, always adding `--nopreopt` when constraints are active and ensuring validation of optional flags.
   - Validate implicit solvent settings against the supported combinations table below before appending `--alpb` or `--gbsa` flags, raising informative errors otherwise.
   - Support custom scratch directory management and allow user override for crest binary path.

3. **Execution Layer**
   - Use Python’s `subprocess` to launch CREST, capturing stdout/stderr for logging.
   - Stream logs to an optional callback or progress printer; raise descriptive errors on non-zero exit codes.
   - Clean up temporary files unless `keep_temps=True` is specified.

4. **Result Handling**
   - Collect output files (`crest_best.xyz`, `full_ensemble.xyz`, populations) and load them back into ASE `Atoms` or trajectory objects.
   - Provide helper functions to parse ensembles and filter conformers by energy or population thresholds.
   - Return a structured result object that surfaces metadata (wall time, configuration, paths).

5. **User-Facing Utilities**
   - Add convenience wrappers for common workflows (single-run microsolvation, constrained runs derived from atom index lists or automatic solute-wide constraints).
   - Provide query helpers (e.g., `list_supported_implicit_solvents`, `supports_implicit_solvent`) backed by the table below.
   - Integrate minimal logging configuration compatible with user-defined loggers.

6. **Testing & Validation**
   - Implement unit tests with fixtures for mock CREST executions and sample output parsing.
   - Include integration tests gated behind an environment variable that requires CREST to be installed.
   - Ensure compatibility with Python 3.9+ and ASE conventions (e.g., `ase.io` readers/writers).

## Usage Examples

### Python API Workflow

```python
from ase.io import read
from microsolvator import Microsolvator, MicrosolvatorConfig

solute = read("benzoic_acid.xyz")
solvent = read("water.xyz")

config = MicrosolvatorConfig(
    nsolv=3,
    temperature=298.0,
    mdtime=50.0,
    crest_executable="crest",
    implicit_solvent="water",
    nopreopt=False,
)

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    keep_temps=False,
)

lowest_cluster = result.best_structure  # ASE Atoms object
ensemble = result.ensemble  # list[Atoms]
```

### Supplying Constraints via Atom Indices

```python
from microsolvator import Microsolvator

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrained_indices=[3, 4, 6, 7, 8, 16, 17, 18, 19, 20, 21, 22],
)
```

`Microsolvator` will emit an `.xcontrol` file for the provided indices, append `--nopreopt` automatically, and clean the temporary file unless `keep_temps=True`.

### Constraining the Entire Solute

```python
result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
    constrain_solute=True,
)
```

When `constrain_solute=True`, the wrapper constrains every atom in the solute structure and activates `--nopreopt` to preserve the input geometry.


### Checking Implicit Solvent Support

```python
from microsolvator import Microsolvator

if not Microsolvator.supports_implicit_solvent(method="gfn2", model="alpb", solvent="water"):
    raise ValueError("Requested implicit solvent combination is unavailable")

result = Microsolvator.run(
    solute=solute,
    solvent=solvent,
    config=config,
)
```

Helper methods return booleans or structured listings sourced from the table below so that invalid combinations are caught before launching CREST.


### Implicit Solvents

```python
solvent_data = {
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
            "h2o": True
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
            "h2o": True
        }
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
            "h2o": True
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
            "h2o": True
        }
    },
    "gfn-ff": {
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
        "h2o": True
    }
}
```