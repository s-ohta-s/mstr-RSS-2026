# Reproduction Materials: Multivariate Spatio-Temporal Regression with Penalized Model Selection and an Empirical Application

Code, data, and outputs for the synthetic-data experiments and the Kansai 198-municipality application of the multivariate general nesting spatio-temporal (MGNST) model and its nested submodels (MSAR, MSEM, and eight simpler candidates), with penalized information criteria pAIC/pBIC for model selection.

> **Paper:** Ryuei Nishii, Saeko Ohta and Shojiro Tanaka (2026)  Multivariate Spatio-Temporal Regression with Penalized Model Selection and an Empirical Application

All result files referenced in the paper are included in the `output/` directories, so the results can be inspected without re-running anything. Re-running is possible with the commands below but is computationally heavy (the n = 400 synthetic-data run took ~21 hours on 8 workers).

## Repository layout

```
s5-0-Three_Implementations/ three separately developed implementations (A/B/C)
                              of the estimation algorithm, with their outputs;
                              see MSTR_implementations_README.md
  estimation-algorithm-mstr.pdf   estimation algorithm document
  simulated_data_1111_n400_T5.csv   shared input panel (n = 400, K = 2, T = 5)
  spatial_weights_n400.csv        shared spatial weight matrix
  implementation_A/implement-A.r  reference implementation; output_n400/
  implementation_B/implement-B.r  instrumented unified runner; results_n400/
  implementation_C/implement-C.r  compact direct implementation; xlsx/csv outputs
  calc_environment/               benchmark platform and R version records
s5-1-Synthetic_Data_Estimation/  synthetic-data experiments
  scripts/run_mc.R            experiment driver (Experiments A and B)
  scripts/plot_frobenius_boxplots_v2.R   figure generation
  R/                          estimation and simulation functions
  data/simulated/             spatial weight matrices (n = 100, 400, 900)
  output/nXXX/                results per n; see output/column_descriptions.md
  output/figures_v2/          figures used in the paper
s5-2-Kansai_Estimation/     real-data application (198 Kansai municipalities)
  scripts/run_real_data.r     estimation, γ search, comparison tables
  scripts/run_identifiability_real.r   identifiability diagnostics
  R/                          estimation functions (diverged from
                              s5-1-Synthetic_Data_Estimation/R)
  data/real/                  input data; see data_source_and_processing_en.md
  output/                     tables and residuals; see column_descriptions.md
s5-2-Kansai_Map/            choropleth maps of fitted values and residuals
  create_heatmaps_en.r        map script (shapefiles not included; see below)
  region_id_mapping198_en.csv region ID to municipality name mapping
  heatmaps/                   generated maps and scatter plots (PNG)
```

**Status of `s5-0-Three_Implementations/`.** This directory is supporting
reference material for the paper, not a part of it: none of its contents are
cited in the paper, and none of the results reported in the paper are produced
by it. The three implementations (A/B/C) were developed separately from each
other and are provided so that readers can inspect how the estimation algorithm
of the paper can be realized in practice, and see the cross-check of their
numerical agreement documented in `MSTR_implementations_README.md`. It uses its
own bundled input files and is therefore independent of the directories below;
it is not part of the reproduction sequence.

## Reproduction order

The four steps below reproduce the results reported in the paper.
`s5-0-Three_Implementations/` is outside this sequence (see the note above).

1. `s5-1-Synthetic_Data_Estimation/scripts/run_mc.R` — synthetic-data experiments
   (`--outdir` selects the output directory; reads the weight matrices in
   `s5-1-Synthetic_Data_Estimation/data/simulated/`)
2. `s5-1-Synthetic_Data_Estimation/scripts/plot_frobenius_boxplots_v2.R` — figures
   (reads each `output/nXXX/results_A.csv`)
3. `s5-2-Kansai_Estimation/scripts/run_real_data.r`, then
   `s5-2-Kansai_Estimation/scripts/run_identifiability_real.r`
   (the latter depends on the former's output)
4. `s5-2-Kansai_Map/create_heatmaps_en.r` — maps
   (depends on the `fitted_residuals_*` output of step 3)

> **Note for step 4.** The input paths at the top of `create_heatmaps_en.r`
> still refer to the previous directory name and must be edited before the
> script is run: `RESID_DIR` (§0, line 16), `PATH_DATA` (line 17) and `PATH_W`
> (line 19) point to `../Kansai_Estimation/...` and should be changed to
> `../s5-2-Kansai_Estimation/...`. `PATH_MAPPING` and `OUTPUT_DIR` are relative
> to the script's own directory and need no change.

## Environment and reproducibility

- R 4.5.2 for the reproduction steps above. The three implementations in
  `s5-0-Three_Implementations/` were run under R 4.4.2; their platform and
  version records are in `s5-0-Three_Implementations/calc_environment/` and
  §6 of `s5-0-Three_Implementations/MSTR_implementations_README.md`.
- Synthetic-data run seeds and settings are recorded in each
  `s5-1-Synthetic_Data_Estimation/output/nXXX/mc_run_info.txt`
  (base_seed = 20260711, 300 replications per cell, 0/3600 failures at n = 400)
- Packages: see the `library()`/`requireNamespace()` calls at the top of each  script. The map script additionally requires
  sf, ggplot2, dplyr, tidyr, viridis, patchwork, stringr, spdep, rmapshaper.

## Data

- Real data: `s5-2-Kansai_Estimation/data/real/` — sources and processing are  documented in `data_source_and_processing_en.md` in that directory.
- Shapefiles for the maps are **not included**. Download URLs
  (MLIT National Land Numerical Information: N03 administrative boundaries,  W09 lakes) are given in §1 of `s5-2-Kansai_Map/create_heatmaps_en.r`.

## License

MIT (see `LICENSE`). Note that the license covers the code and documents in this repository; the underlying statistical data and map data remain subject to the terms of their original providers (see `s5-2-Kansai_Estimation/data/real/data_source_and_processing_en.md`).
