# Beyond straight lines: migration costs considering geography enhance tracing human genetic ancestry

This repo contains the tex files (under ms) as well as SLiM and R files (under data) used for the simulations, analyses, and visualizations associated with the [manuscript](https) describing the gaia software.

The full workflow is in `code/master.R`, which sources every other script and reproduces the analysis from raw inputs to final figures and tables.

## Overview

The project incldes:

- A heuristic 1d example of tree sequence structure, demonstrating why MPR algorithm works in general (figure 1).
- Comparison of the practicality of applying heterogeneous migration costs across multiple network dimensions (figure 2).
- A simulated case study of human migrations out of Africa, comparing the proposed friction method with the naive baseline (figure 3).
- A real world example of human migrations into the Americas, verifying if the estimation results avoid major geographical barriers (figure 4).

The code produces:

- intermediate results of tables and trees,
- summary tables,
- manuscript figures,
- supplementary figures.

## Repository structure

```text
geographical_ancestry_inference-paper/
├── code/
│   ├── master.R
│   ├── functions.R
│   ├── generation/
│   ├── simulation/
│   ├── visualization/
├── data/
│   ├── flux/
│   ├── genetics/
│   │   ├── subsets/
│   ├── geo/
│   ├── math/
│   ├── mpr/
│   ├── pop/
│   ├── trees/
├── output/
│   ├── figures/
│   ├── tables/
└── README.md
```

## How to reproduce the analysis

1. Open the repository in RStudio or VS Code.
2. Set the working directory to the project root.
3. Load the following input files:

- in data/geo: 2020_walking_only_friction_surface.geotiff obtained from [the Malaria Atlas Project](https://malariaatlas.org/project-resources/accessibility-to-healthcare/); landgrid_adjmat_naive_afro-eurasia.csv and landgrid_afro-eurasia.gpkg obtained from [GAIA paper](https://github.com/blueraleigh/gaia)
- in data/pop: popc_5000BC.asc (and other files including 10000BC, 0AD, 1500AD, 2000AD) obtained from [History Database of the Global Environment (HYDE)](https://geo.public.data.uu.nl/vault-hyde/hyde35_c9_apr2025%5B1749214444%5D/original/gbc2025_7apr_base/zip/)
- in data/genetics: hgdp_tgp_sgdp_high_cov_ancients_chr18_p.dated.trees obtained from [Unified Genealogy of Modern and Ancient Human Genomes dataset](https://zenodo.org/records/5512994)

4. Install SLiM software from https://messerlab.org/slim/
5. Run:

```r
source("code/master.R")
```

`master.R` will:

- load required packages,
- create the output folders if needed,
- run data preparation,
- run estimation of ancestral spatial locations via GAIA,
- compare the friction and naive method,
- and generate the figures and tables.

## Citation

If you use this code, please cite the associated manuscript and acknowledge the repository.
