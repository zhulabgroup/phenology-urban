# phenology-urban

Urban tree phenology workflow and analysis code.

## Data source

Project data are hosted on Zenodo:

- https://doi.org/10.5281/zenodo.21709425

After downloading, place the data under the `data/` folder (matching this repository structure).

## Workflow location

The full workflow description is documented in:

- `vignettes/outline.qmd`

This document summarizes data processing, modeling, and figure generation for the project. A rendered version is also available [here](https://public.websites.umich.edu/~jializ/urban-phe/outline.html).

## Pipeline overview

Main scripts are organized by step under `scripts/`:

1. `scripts/01_tidy_tree_inventory/`
	- Collect, clean, and summarize urban tree inventories.
	- Assign leaf phenology type from TRY data.

2. `scripts/02_process_climate/`
	- Process UHAT, Daymet climate data.
	- Build tree-level climate inputs.

3. `scripts/03_process_PS/`
	- Download/process PlanetScope imagery.
	- Extract individual-tree phenology metrics (SOS/EOS).

4. `scripts/04_stats_model/`
	- Fit variance partitioning, spatial models, and PGLMM models.

5. `scripts/05_plotting/`
	- Generate manuscript figures.

## Key data folders

- `data/input/`: processed inputs (tree list, climate, phenology)
- `data/intermediate/`: model intermediate outputs
- `data/output/`: final model outputs and summaries

## Notes

- This repository assumes the data have been downloaded from Zenodo first.
- See `vignettes/outline.qmd` for detailed code snippets and figure-specific workflows.