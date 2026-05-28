# geographical_ancestry_inference-paper

This repo contains the tex files (under ms) as well as SLiM and R files (under data) used for the simulations, analyses, and visualizations associated with the [manuscript](https) describing the gaia software.

The full workflow is in `code/master.R`, which sources every other script and reproduces the analysis from raw inputs to final figures and tables.

## Overview

The project examines:

- where Chinese peacekeeping deployments are located relative to other troop contributing countries,
- how close deployments are to **primary roads** and **large cities**,
- how spatial coverage and overlap evolve over time,
- and how these patterns compare with other contributors, especially **P4** and **non-P5** countries.

The code produces:

- raster outputs for spatial coverage and overlap,
- summary tables,
- manuscript-ready figures,
- supplementary figures,
- and a case-study map for South Sudan.

## Repository structure

```text
geographical_ancestry_inference-paper/
├── code/
│   ├── master.R
│   ├── generation/
│   │   ├── create_landgrid.R
│   │   ├── 
│   ├── simulation/
│   │   ├── 
│   │   ├── 
│   ├── visualization/
│   │   ├── visualize_line.R
│   │   ├── tree_demo.R
├── data/
├── output/
└── README.md
```

## How to reproduce the analysis

1. Open the repository in RStudio or VS Code.
2. Set the working directory to the project root.
3. Run:

```r
source("code/master.R")
```

`master.R` will:

- load required packages,
- create the output folders if needed,
- run data preparation,
- compute spatial metrics,
- and generate the figures and tables.

## Main pipeline

### `code/master.R`

Master script that runs the full workflow in order:

1. load packages,
2. create output folders,
3. source helper functions,
4. prepare the data,
5. compute travel-time and spatial-coverage metrics,
6. generate all figures.

### `code/functions.R`

Core helper functions used throughout the project. These include:

- data cleaning and reshaping utilities,
- UN troop aggregation functions,
- mission recoding helpers,
- parallel travel-time routines,
- raster plotting helpers,
- summary/annotation plotting functions used in the manuscript figures.

### `code/dataprep.R`

Builds the analysis-ready dataset. In particular, it:

- loads and filters the UN peacekeeping deployment data,
- keeps cases where China deployed peacekeepers,
- splits the data into China, P5, and non-P5 groups,
- standardizes country codes,
- and saves the cleaned data to `output/rdata/pkopt_by_country.Rdata`.

### `code/traveltime.R`

Computes the core spatial metrics:

- travel time from deployments to spatial objects,
- spatial coverage,
- overlap ratio,
- overlap intensity.

It saves:

- `output/rdata/traveltime_output.Rdata`
- `output/tables/coverage_metrics.csv`
- and additional raster outputs used in the maps.

### `code/coverageplot.R`

Creates the time-series figures for:

- spatial coverage,
- overlap ratio,
- overlap intensity.

It uses summary metrics aggregated by year, TCC group, and sample.

### `code/time2objectsplot.R`

Creates the comparison plots for travel time from deployments to:

- the nearest **primary road**,
- the nearest **large city**.

It also generates yearly trend plots for Chinese deployments and for the pooled contributor groups.

### `code/mappko.R`

Generates the multi-panel map of peacekeeping deployments and spatial coverage constraints for countries where China deployed peacekeepers. It also exports supporting CSV summaries for the figure.

### `code/casestudyplot.R`

Builds the South Sudan case study figure, showing:

- spatial overlap,
- proximity to roads,
- cities above the population threshold,
- and strategic locations such as the capital, embassy, and airport.

### `code/waffletroop_plot.R`

Produces waffle plots summarizing:

- troop size across troop contributing countries,
- fatalities associated with UN peace operations.

### `code/waffle_nonciv_civ_plot.R`

Creates the monthly line/tile-style figure for Chinese peacekeeping personnel composition across missions.

## Inputs

The scripts expect the `input/` folder to contain the raw and auxiliary files used in the analysis, including:

- UN peacekeeping deployment data,
- friction/travel-cost rasters,
- road data,
- populated places,
- mission datasets,
- and any other external spatial layers referenced by the scripts.

Some inputs are downloaded or created during the workflow when missing, while others must already be present in `input/`.

## Outputs

The workflow writes results to `output/`, including:

- `output/rdata/` — intermediate and final `.Rdata` files,
- `output/tables/` — CSV summary tables,
- `output/raster/` — raster outputs for coverage and overlap,
- `output/figsi/` — intermediate/supplementary figures,
- `output/figsmanuscript/` — manuscript-ready figures.

## Notes

- The code assumes it is run from the repository root.
- Spatial operations use `sf`, `terra`, `raster`, and `gdistance`.
- Several scripts use parallel processing, so runtime depends on available cores and data size.
- `sf_use_s2(FALSE)` is used in multiple scripts to avoid topology issues in some polygon operations.

## Citation

If you use this code, please cite the associated manuscript and acknowledge the repository.
