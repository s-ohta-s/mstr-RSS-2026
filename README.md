# Reproduction Materials: Multivariate Spatio-Temporal Regression with Penalized Model Selection and an Empirical Application

Code, data, and outputs for the Monte Carlo experiments and the Kansai 198-municipality application of the multivariate general nesting spatio-temporal (MGNST) model and its nested submodels (MSAR, MSEM, and eight simpler candidates), with penalized information criteria pAIC/pBIC for model selection.

> **Paper:** Ryuei Nishii, Saeko Ohta and Shojiro Tanaka (2026)  Multivariate Spatio-Temporal Regression with Penalized Model Selection and an Empirical Application

All result files referenced in the paper are included in the `output/` directories, so the results can be inspected without re-running anything. Re-running is possible with the commands below but is computationally heavy (the n = 400 Monte Carlo run took ~21 hours on 8 workers).

## Repository layout

```
montecarlo/                 Monte Carlo experiments
  scripts/run_mc.R            experiment driver (Experiments A and B)
  scripts/plot_frobenius_boxplots_v2.R   figure generation
  R/                          estimation and simulation functions
  data/simulated/             spatial weight matrices (n = 100, 400, 900)
  output/nXXX/                results per n; see output/column_descriptions.md
  output/figures_v2/          figures used in the paper
Kansai_Estimation/          real-data application (198 Kansai municipalities)
  scripts/run_real_data.r     estimation, γ search, comparison tables
  scripts/run_identifiability_real.r   identifiability diagnostics
  R/                          estimation functions (diverged from montecarlo/R)
  data/real/                  input data; see data_source_and_processing_en.md
  output/                     tables and residuals; see column_descriptions.md
Kansai_Map/                 choropleth maps of fitted values and residuals
  create_heatmaps_en.r        map script (shapefiles not included; see below)
Three_Reference_Implementations/   three separately developed verification
  implementations (A/B/C) with their outputs; see README_implementations.md
```

## Reproduction order

1. `montecarlo/scripts/run_mc.R` — Monte Carlo experiments
   (`--outdir` selects the output directory; reads the weight matrices in
   `montecarlo/data/simulated/`)
2. `montecarlo/scripts/plot_frobenius_boxplots_v2.R` — figures
   (reads each `output/nXXX/results_A.csv`)
3. `Kansai_Estimation/scripts/run_real_data.r`, then
   `Kansai_Estimation/scripts/run_identifiability_real.r`
   (the latter depends on the former's output)
4. `Kansai_Map/create_heatmaps_en.r` — maps
   (depends on the `fitted_residuals_*` output of step 3)

## Environment and reproducibility

- R 4.5.2 (versions used for the reference implementations are documented in
  `Three_Reference_Implementations/calc_environment/`)
- Monte Carlo seeds and settings are recorded in each
  `montecarlo/output/nXXX/mc_run_info.txt`
  (base_seed = 20260711, 300 replications per cell, 0/3600 failures at n = 400)
- Packages: see the `library()`/`requireNamespace()` calls at the top of each  script. The map script additionally requires
  sf, ggplot2, dplyr, tidyr, viridis, patchwork, stringr, spdep, rmapshaper.

## Data

- Real data: `Kansai_Estimation/data/real/` — sources and processing are  documented in `data_source_and_processing_en.md` in that directory.
- Shapefiles for the maps are **not included**. Download URLs
  (MLIT National Land Numerical Information: N03 administrative boundaries,  W09 lakes) are given in §1 of `Kansai_Map/create_heatmaps_en.r`.

## License

MIT (see `LICENSE`). Note that the license covers the code and documents in this repository; the underlying statistical data and map data remain subject to the terms of their original providers (see `Kansai_Estimation/data/real/data_source_and_processing_en.md`).
