# phenology-urban

Urban tree phenology workflow and analysis code.

## Data source

Project data are hosted on Zenodo:

- https://doi.org/10.5281/zenodo.21709425

After downloading, place the data under the `release/` folder (matching this repository structure).

## Workflow location

The full workflow description is documented in:

- `doc/outline.qmd`

This document summarizes data processing, modeling, and figure generation for the project.

## Pipeline overview

Main scripts are organized by step under `script/`:

1. `script/01_tidy_tree_inventory/`
	- Collect, clean, and summarize urban tree inventories.
	- Assign leaf phenology type from TRY data.

2. `script/02_process_climate/`
	- Process UHAT, Daymet climate data.
	- Build tree-level climate inputs.

3. `script/03_process_PS/`
	- Download/process PlanetScope imagery.
	- Extract individual-tree phenology metrics (SOS/EOS).

4. `script/04_stats_model/`
	- Fit variance partitioning, spatial models, and PGLMM models.

5. `script/05_plotting/`
	- Generate manuscript figures.

## Key data folders

- `release/input/`: processed inputs (tree list, climate, phenology)
- `release/intermediate/`: model intermediate outputs
- `release/output/`: final model outputs and summaries

## Notes

- This repository assumes the data have been downloaded from Zenodo first.
- See `doc/outline.qmd` for detailed code snippets and figure-specific workflows.