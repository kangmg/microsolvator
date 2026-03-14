# Installation

## Package

```bash
pip install git+https://github.com/kangmg/microsolvator.git
```

## External binaries

`microsolvator` requires **CREST** and **xTB**. Use the built-in installers or provide your own paths.

=== "Built-in installer"

    ```python
    from microsolvator import install_crest, install_xtb

    install_crest()
    install_xtb()
    ```

    Binaries are stored in `microsolvator/_bin/` and resolved automatically.

=== "Manual paths"

    ```python
    from microsolvator import MicrosolvatorConfig

    config = MicrosolvatorConfig(
        nsolv=5,
        crest_executable="/opt/crest/crest",
        xtb_executable="/opt/xtb/bin/xtb",
    )
    ```

=== "Environment variables"

    ```bash
    export CREST_BIN=/opt/crest/crest
    export XTB_BIN=/opt/xtb/bin/xtb
    ```

## Binary resolution order

For both CREST and xTB, the resolver checks in this order:

1. Explicit path in `MicrosolvatorConfig`
2. `CREST_BIN` / `XTB_BIN` environment variables
3. `microsolvator/_bin/` (installed via `install_crest` / `install_xtb`)
4. System `PATH`

## Workflow dependencies

The `microsolvator.workflow` subpackage additionally requires:

- **Packmol** — for box solvation (Step 2)
- An **ASE calculator** — for MD equilibration and relaxation (Steps 3–4)

```bash
# Packmol (example: conda)
conda install -c conda-forge packmol

# Calculator example: xtb-python
pip install xtb-python
```

Packmol is resolved from `PackmolConfig.packmol_executable` or system `PATH`.
