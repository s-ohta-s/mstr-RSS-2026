# Column Descriptions for the Monte Carlo Output Files

This note documents the columns of `results_A.csv`, `results_B.csv`, and
`identifiability_diag.csv` found in each `nXXX/` subdirectory (n = 100, 400, 900).
Each `nXXX/` directory also contains `mc_run_info.txt` (run configuration and
seeds), `mc_timings.csv` (per-replication timing), and
`true_parameters_<dgp>_<paramset>.csv` (true parameter values used by the DGP).

## Model ID convention

Model IDs are four-character codes giving the status of the blocks
(**R**, **Λ**, **A**, **Σ**) — spatial lag, spatial error, temporal lag,
error covariance:

- `1` — block unrestricted (fully estimated)
- `0` — block absent (fixed to zero)
- `d` — block restricted to diagonal (independent univariate equations)

The six γ-penalized "full" models are `1111` (MGNST), `0111` (MSEM),
`1011` (MSAR) with temporal lag, and `1101`, `0101`, `1001` (the same three
without temporal lag). The five simple candidates are `0011` (VARX),
`000d` (independent regression), `d0dd` (independent SAR), `0ddd`
(independent SEM), and `dddd` (independent GNS).

## `results_A.csv` — Experiment A (estimation accuracy under the true model)

One row per parameter per replication. The true model is fitted with a
single-stage γ search; estimates and SEs refer to the fit at the
**pAIC-selected γ** unless noted otherwise.

| Column | Description |
|---|---|
| `dgp` | True data-generating process, also the fitted model (`1111`, `0111`, or `1011`) |
| `paramset` | True-parameter scenario: `weak`, `mid`, `strong` (the paper's Weak / Moderate / Strong levels), or `realistic` — an additional scenario based on real-data estimates that is not reported in the paper (the paper's 8,100 replications correspond to the weak/mid/strong cells only) |
| `n` | Number of regions |
| `rep` | Replication number (1–300) |
| `seed` | RNG seed used to generate this replication's dataset |
| `parameter` | Parameter name: `beta_intercept_*`, `beta_<x>_<y>`, `A[i,j]` (temporal lag), `R[i,j]`, `Lambda[i,j]`, `Sigma[i,j]` |
| `true` | True value (see `true_parameters_<dgp>_<paramset>.csv`) |
| `estimate` | Estimate at the pAIC-selected γ |
| `se_psi` | Standard error used for the coverage indicators below: Ψ-based (paper Eq. (18), GLS covariance) for the β/A rows, penalized-profile-Hessian SE for the R/Λ rows |
| `se_hessian` | Full-Hessian sandwich SE; computed only when the expensive full-Hessian option is enabled — `NA` in these runs |
| `covered95_psi` | 1 if the true value lies in the 95% normal CI built from `se_psi`, 0 otherwise, `NA` if the SE is unavailable |
| `covered95_hess` | Same using `se_hessian` (`NA` in these runs) |
| `covered95_main` | Main coverage indicator of this experiment (identical to `covered95_psi`); per-parameter coverage is not itself tabulated in the paper |
| `gamma` | γ selected by pAIC |
| `gamma_pBIC` | γ selected by pBIC |
| `estimate_pBIC` | Estimate at the pBIC-selected γ (point estimate only; no SE) |
| `df_eff` | Effective degrees of freedom at the selected γ (paper Eq. (20)) |
| `loglik` | Log-likelihood of the pAIC-selected fit |
| `converged` | Optimizer convergence flag |
| `stationarity_ok` | Overall validity of the estimate: spatial invertibility (`min_det_*` > 0), positive-definite Σ̂, and — where applicable — temporal stability (`rho_dyn` < 1) |
| `rho_dyn` | Spectral radius of the estimated space–time dynamics; < 1 indicates a stationary process |
| `min_det_R` | min over eigenvalues ω of W of det(I − ω R̂); > 0 certifies invertibility of (I − R̂⊗W) |
| `min_det_Lambda` | Same for Λ̂ |
| `prop1_rank_ok` | Whether the Proposition 1 rank condition (identifiability of R and β) holds on the realized design, evaluated at the true parameters |
| `prop1_min_sv` | Margin of that condition: minimum singular value of the column-normalized matrix [G_µ0, X] |
| `hess_min_eig_gamma0` | Minimum eigenvalue of the unpenalized profile Hessian at γ = 0 (local-identification diagnostic) |

Columns from `gamma` onward are replication-level quantities and repeat on
every row belonging to the same replication.

## `results_B.csv` — Experiment B (model selection among 11 candidates)

One row per replication. All 11 candidate models are fitted; the penalized
criteria pAIC/pBIC of each full model are minimized over the γ grid, and the
simple models are evaluated at γ = 0 (so their pAIC/pBIC equal ordinary
AIC/BIC).

| Column | Description |
|---|---|
| `dgp`, `paramset`, `n`, `rep`, `seed` | Cell identifiers, as in `results_A.csv` |
| `selected_pAIC` | Model ID with the smallest pAIC |
| `selected_pBIC` | Model ID with the smallest pBIC |
| `true_model` | True DGP (equals `dgp`) |
| `hit_pAIC` | 1 if `selected_pAIC` equals the true model, else 0 |
| `hit_pBIC` | 1 if `selected_pBIC` equals the true model, else 0 |
| `n_models_ok` | Number of the 11 candidates that were fitted successfully |
| `pAIC_<id>`, `pBIC_<id>` | Criterion values for candidate `<id>` (`NA` if that fit failed) |

## `identifiability_diag.csv` — identifiability diagnostics

Numerical checks of the two sufficient identifiability conditions of the
paper (§3.2). Two kinds of rows, distinguished by `scope`:

- `scope = "W"` — one row per weight matrix (per n), computed once.
  Only `prop2_min_eig` and `W_symmetric` are filled.
- `scope = "cell"` — one row per (dgp, paramset, n) cell, aggregating the
  per-replication Proposition 1 diagnostics from `results_A.csv`.
  Only the `prop1_*` columns are filled.

| Column | Description |
|---|---|
| `scope` | `"W"` or `"cell"` (see above) |
| `n` | Number of regions |
| `dgp`, `paramset` | Cell identifiers (`NA` on `W` rows) |
| `prop1_pass_rate` | Share of replications in which the Proposition 1 rank condition held |
| `prop1_min_sv_min` | Minimum over replications of `prop1_min_sv` |
| `prop1_min_sv_median` | Median over replications of `prop1_min_sv` |
| `prop2_min_eig` | Minimum eigenvalue of the normalized Gram matrix of {I, W, W′, W′W}; > 0 certifies their linear independence (Proposition 2). Note: the paper reports the minimum singular value of the normalized [vec(I), vec(W), vec(W′), vec(W′W)] matrix instead, which is the square root of this eigenvalue |
| `W_symmetric` | Whether the row-standardized W is symmetric (tolerance 1e-12); the Proposition 2 condition concerns the generic nonsymmetric case |
