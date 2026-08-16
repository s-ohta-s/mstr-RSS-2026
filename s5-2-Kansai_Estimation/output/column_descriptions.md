# Column Descriptions for the Real-Data Comparison Tables

This note documents `comparison_table_pAIC.csv` and `comparison_table_pBIC.csv`,
the cross-model parameter comparison tables for the Kansai 198-municipality
application. The two files have identical structure; they differ only in the
criterion used to select γ for each of the six penalized full models
(pAIC vs. pBIC). Both are written by `scripts/run_real_data.r`.

## Layout

One row per parameter (or fit statistic), one **group of five columns per
candidate model**:

```
estimate.<id>, se_psi.<id>, se_hessian.<id>, signif_psi.<id>, signif_hessian.<id>
```

where `<id>` is the model ID. Model IDs are four-character codes giving the
status of the blocks (**R**, **Λ**, **A**, **Σ**) — spatial lag, spatial
error, temporal lag, error covariance — with `1` = unrestricted, `0` = absent,
`d` = diagonal (independent univariate equations):

- Full (γ-penalized) models: `1111` (MGNST), `0111` (MSEM), `1011` (MSAR),
  and the no-temporal-lag versions `1101`, `0101`, `1001`
- Simple models (γ = 0): `0011` (VARX), `000d` (independent regression),
  `d0dd` (independent SAR), `0ddd` (independent SEM), `dddd` (independent GNS)

An empty cell means the quantity is not defined for that model (e.g. the
parameter is not part of the model, or that type of SE is not computed for it).

## Rows

| Row | Description |
|---|---|
| `beta_intercept_<y>`, `beta_<x>_<y>` | Regression coefficients for response `<y>` (`y1`, `y2`) |
| `A[i,j]` | Temporal-lag (AR) coefficients |
| `R[i,j]` | Spatial-lag coefficient matrix |
| `Lambda[i,j]` | Spatial-error coefficient matrix |
| `Sigma[i,j]` | Error covariance (upper triangle only) |
| `loglik`, `AIC`, `BIC` | Log-likelihood and ordinary information criteria |
| `pAIC`, `pBIC` | Penalized information criteria at the selected γ (equal AIC/BIC for the simple models) |
| `d_eff` | Effective degrees of freedom at the selected γ (paper Eq. (20)) |
| `gamma` | Selected ridge parameter γ (0 for the simple models) |
| `n_params` | Nominal number of estimated parameters |
| `pseudo_R2`, `adj_R2` | Pseudo-R² = corr(y, ŷ)² using the trend prediction (MSAR/MGNS: ŷ = (I − R̂⊗W)⁻¹Xβ̂; MSEM/VARX/OLS: ŷ = Xβ̂), and its adjusted version |
| `pseudo_R2_y1`, `pseudo_R2_y2` | Per-response pseudo-R² |

Fit-statistic rows carry values only in the `estimate.<id>` column of each group.

## The five columns per model

| Column | Description |
|---|---|
| `estimate.<id>` | Point estimate at the selected γ |
| `se_psi.<id>` | Ψ-based (GLS covariance) standard error — defined for the β/A rows of the full models; for the diagonal models it is the per-equation SE reported by `spatialreg`. Blank on R/Λ rows |
| `se_hessian.<id>` | Full-Hessian standard error (available for the six full models; covers β/A and R/Λ rows). Blank for the simple models |
| `signif_psi.<id>` | Significance code from the two-sided normal test estimate/`se_psi` |
| `signif_hessian.<id>` | Same using `se_hessian` |

Significance codes: `***` p < 0.001, `**` p < 0.01, `*` p < 0.05,
`.` p < 0.1, blank otherwise (or when the SE is unavailable).

Relation to Table 3 of the paper: for regression coefficients the paper
assesses significance conservatively with the larger of the two standard
errors — i.e. its symbol corresponds to the weaker of `signif_psi` and
`signif_hessian` — while for the spatial parameters (R/Λ rows) it uses the
Hessian-based SE (`signif_hessian`).
