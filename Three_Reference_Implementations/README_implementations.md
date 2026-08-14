# Three Reference Implementations of Simultaneous Parameter Estimation for Multivariate Spatio-Temporal Regression (MSTR)

This directory provides three separately developed, related R implementations — `implement-A.r`, `implement-B.r`, and `implement-C.r` — of the maximum likelihood estimation procedure for the multivariate general nesting spatio-temporal (MGNST) model and its ten nested submodels, as formulated in the paper accompanying this repository (see the top-level README for the citation; all equation and section numbers below refer to that paper: model Eqs. (4)–(7), likelihood and GLS estimation Eqs. (11)–(18), penalized estimation and effective degrees of freedom Eqs. (19)–(21) and §4.4, log-determinants Appendix B, identifiability §3.2 and Appendix A). All three follow the four-stage scheme of §4.2 of the paper: (S1) univariate-GNS initial values, (S2) inner iterative GLS updates of $\hat{\boldsymbol{\beta}}$ (Eq. 15) and $\hat{\boldsymbol{\Sigma}}$ (Eq. 17), (S3) the profile log-likelihood in the spatial parameters (step S3 of §4.2), and (S4) quasi-Newton updates of the free entries of $(\boldsymbol{R}, \boldsymbol{\Lambda})$, with penalized information criteria pAIC/pBIC (§4.4) based on the effective degrees of freedom of Eq. (20).

The three programs are not fully independent — they implement the same specification and share design conventions — but they were written separately and differ in code base, log-determinant scheme, optimization strategy, and operational design. A provides a partially independent numerical reference (general-$K$ symmetric-function log-determinant, external `spatialreg` initialization), while B and C are two architecturally distinct realizations of the $K = 2$ estimator. Agreement among their outputs therefore constitutes a meaningful, though not fully independent, cross-check.

**Scope of the numerical verification.** All verification and benchmark results below concern the conditional estimation at a single target period $t = 2$, with the $t = 1$ observations entering only through the AR(1) regressors of Eq. (7). They do not verify a joint likelihood over all five periods of the supplied panel. (The dataset contains $T = 5$ periods; `TARGET_TIME`/`time_now` is set to 2 in all three scripts.)

## 1. The three implementations

### `implement-A.r` — Reference implementation (general $K$)

A comprehensive development version (13,994 lines) integrating all modules produced during development: data preparation, eleven model-specific fitters, penalized spatial estimation, CSV output, and a parallel experiment runner. Its distinguishing feature is a faithful implementation of the log-determinant $\log\lvert\boldsymbol{I} - \boldsymbol{R}\otimes\boldsymbol{W}\rvert$ via the elementary-symmetric-polynomial expansion of Appendix B (Eqs. (27)–(28)), computed from power-sum traces $r_j = \operatorname{tr}(\boldsymbol{R}^j)$ by Newton's identities and valid for general response dimension $K \leq 7$. Admissibility of candidate $(\boldsymbol{R}, \boldsymbol{\Lambda})$ is checked explicitly by a `spectral_radius()` routine that tests $\min_i\det(\boldsymbol{I}_K - \omega_i\boldsymbol{M}) > 0$ over all eigenvalues $\omega_i$ of $\boldsymbol{W}$ and reports a graded violation measure. Initial values are obtained through `spatialreg::lagsarlm` / `errorsarlm`, so A carries the largest dependency footprint (`spdep`, `spatialreg`, `Matrix`, `numDeriv`, `parallel`). The $\gamma$ search is two-stage (coarse grid, then log-scale refinement around the coarse optimum). A is retained as (i) the reference for extending the estimator beyond $K = 2$ and (ii) a partially independent implementation for cross-checking B and C.

### `implement-B.r` — Fully instrumented unified runner ($K = 2$)

A self-contained runner (2,330 lines) specialized to $K = 2$. All eleven model IDs are decoded into masks on $(\boldsymbol{R}, \boldsymbol{\Lambda}, \boldsymbol{A}, \boldsymbol{\Sigma})$ and estimated through a single shared pipeline. Optimization is by multi-start BFGS (`optim`, `reltol = 1e-8`, `maxit = 220` with a 320-iteration retry) on elementwise tanh-transformed variables, with a deterministic multi-start set (univariate initial values, zero starts, sign and scale perturbations) augmented — not replaced — by warm starts along the $\gamma$ path. The log-determinant is evaluated eigenvalue-wise as $\sum_i \log\lvert\det(\boldsymbol{I}_2 - \omega_i\boldsymbol{R})\rvert$, the $K = 2$ form of the determinant identity of Appendix B. Effective degrees of freedom use a numerically evaluated profile Hessian (Section 2.4 below) with a "candidate" screening mode (`DEFF_EVALUATION_MODE = "candidate"`, the default): rigorous lower/upper pAIC/pBIC bounds identify the $\gamma$ grid rows that could still be optimal, and the expensive numerical Hessian is evaluated only for those rows, so the selected pAIC/pBIC values are the same as under exhaustive evaluation. Every fit reports `n_starts` and convergence diagnostics. B also contains the operational hardening previously attributed elsewhere: interactive/environment-variable control of the output root (`MSTR_OUT_ROOT`, `MSTR_USE_CURRENT_DIR`), and safe CSV/workbook writes that fall back to timestamped filenames when files are locked. Dependencies: `Matrix`, `openxlsx`, `parallel` (`numDeriv` optional; `optimHess` fallback).

### `implement-C.r` — Compact direct implementation ($K = 2$)

C (646 lines) is a compact, direct transcription of the estimation scheme. The profile likelihood forms $\boldsymbol{I}_{Kn} - \boldsymbol{R}\otimes\boldsymbol{W}$ and $\boldsymbol{I}_{Kn} - \boldsymbol{\Lambda}\otimes\boldsymbol{W}$ explicitly (dense $800\times800$ for this dataset) for the GLS steps, while the log-determinant is evaluated from the eigenvalues of the $K\times K$ coefficient matrix and of $\boldsymbol{W}$ as $\sum_{i,j}\log\lvert 1 - a_j\omega_i\rvert$ with a strict admissibility check (any real factor $\leq 10^{-12}$ $\Rightarrow$ the candidate is rejected with objective $10^{10}$). Optimization is single-trajectory `optim(method = "L-BFGS-B")` with box constraints $-0.99 \leq \theta \leq 0.99$ on each free spatial parameter and warm starts carried along the $\gamma$ grid (initial value 0.01 for every free parameter at the first grid point; no multi-start). Effective degrees of freedom use a `numDeriv::hessian` of the negative profile log-likelihood at each $\gamma > 0$, with a documented fallback to the nominal spatial parameter count when the Hessian computation fails. Models are parallelized with `doParallel`/`foreach` (one model per worker), and results are written to xlsx. Dependencies: `Matrix`, `parallel`, `doParallel`, `foreach`, `openxlsx`, `numDeriv`.

**Correction relative to an earlier draft of this document.** An earlier draft stated that B and C share a byte-identical numerical core differing only in I/O, and attributed candidate-$\gamma$ screening to C. Inspection of the published files shows the opposite: candidate screening is implemented in B and absent from C, and B and C are architecturally distinct programs (multi-start tanh-space BFGS with screening in B; single-start box-constrained L-BFGS-B with dense Kronecker algebra in C). All statements below are grounded in the published source files; their SHA-256 checksums are listed in Section 6.

## 2. Design comparison

| Criterion | A | B | C |
|---|---|---|---|
| Size / structure | 13,994 lines, all development modules | 2,330 lines, unified runner | 646 lines, compact direct implementation |
| log-det (Appendix B) | Symmetric-polynomial expansion via Newton's identities, general $K \leq 7$ | Eigenvalue-wise $2\times2$ determinants $\sum_i\log\lvert\det(\boldsymbol{I}_2 - \omega_i\boldsymbol{R})\rvert$ | Eigenvalue-pair products $\sum_{i,j}\log\lvert1 - a_j\omega_i\rvert$ with strict positivity check |
| Constraint handling (S4) | Explicit admissibility check $\min_i\det(\boldsymbol{I}_K - \omega_i\boldsymbol{M}) > 0$ | Elementwise tanh ($\lvert\rho_{k\ell}\rvert, \lvert\lambda_{k\ell}\rvert < 1$); singular determinants rejected ($-\infty$) | Box constraints $\pm0.99$; real factors $\leq 10^{-12}$ rejected ($-\infty$) |
| Optimization | `optim` per module (incl. penalized L-BFGS-B); numerical Hessians | Multi-start BFGS in tanh space, warm starts add to start set, retry logic | Single warm-started L-BFGS-B trajectory per $\gamma$ |
| Effective d.o.f. (Eq. (20)) | $\operatorname{tr}[\boldsymbol{I}_{11}(\boldsymbol{I}_{11}+\gamma\boldsymbol{I})^{-1}] + p + K(K+1)/2$, numerical Hessian (GIC module) | Numerical profile $\boldsymbol{\theta}$-Hessian; candidate-$\gamma$ screening with exact-at-selection guarantee; eigenvalues floored at 0 with warning | Numerical profile $\boldsymbol{\theta}$-Hessian per $\gamma > 0$; fallback to nominal count on failure |
| $\gamma$ grid | $\{0\}\cup10^{\{-2,-1.5,\ldots,4\}}$ coarse, then log-scale refinement | $\{0\}\cup10^{\{-2,-1.5,\ldots,4\}}$ (14 points) | $\{0\}\cup10^{\{-2,-1.5,\ldots,4\}}$ (14 points) |
| $\gamma$-penalized models | Full spatial-block models | The six IDs $\{1111, 0111, 1011, 1101, 0101, 1001\}$; others at $\gamma = 0$ with ordinary AIC/BIC | Same six ($\boldsymbol{R}$- or $\boldsymbol{\Lambda}$-code = `"1"`); others at $\gamma = 0$ |
| Initial values (S1) | `spatialreg::lagsarlm` / `errorsarlm` | In-house univariate GNS (cached) | Constant 0.01 per free spatial parameter |
| Parallelism | PSOCK cluster, `MSTR_CORES` (default: physical cores $- 1$), BLAS threads pinned to 1 | Per-model PSOCK; auto worker count (for $n = 400$: $\lfloor0.75(\text{logical cores} - 1)\rfloor$), sequential retry | `doParallel`, $\min(11, \text{logical cores} - 1)$ workers |
| Output | Model-specific CSV writers | xlsx + run-summary / $\gamma$-path CSV, safe writes, output-root prompt/env vars | xlsx summary tables |
| Dependencies | many (`spdep`, `spatialreg`, `systemfit`, …) | few (`Matrix`, `openxlsx`, `parallel`; `numDeriv` optional) | moderate (`Matrix`, `doParallel`, `foreach`, `openxlsx`, `numDeriv`) |

### 2.3 What is actually constrained (stationarity and admissibility)

Step S4 of the paper (§4.2) restricts $(\boldsymbol{R}, \boldsymbol{\Lambda})$ to the admissible region

$$
\mathcal{P} = \{ (\boldsymbol{R}, \boldsymbol{\Lambda}) : \lvert\boldsymbol{I}_{Kn} - \boldsymbol{R}\otimes\boldsymbol{W}\rvert > 0,\ \lvert\boldsymbol{I}_{Kn} - \boldsymbol{\Lambda}\otimes\boldsymbol{W}\rvert > 0 \}.
$$

As implemented, the three programs work with elementwise bounds $\lvert\rho_{k\ell}\rvert < 1$, $\lvert\lambda_{k\ell}\rvert < 1$ — enforced by tanh transformation (B) or L-BFGS-B boxes (C, A's penalized modules) — combined with the admissibility checks described below. Elementwise bounds alone do not guarantee stationarity/invertibility of the matrix process: for non-diagonal $\boldsymbol{R}$, the relevant condition is that all factors $1 - a_j(\boldsymbol{R})\omega_i$ remain bounded away from zero (equivalently, $\det(\boldsymbol{I}_{Kn} - \boldsymbol{R}\otimes\boldsymbol{W}) \neq 0$ along a path connected to $\boldsymbol{R} = \boldsymbol{0}$), and the spectral radius of $\boldsymbol{R}\otimes\boldsymbol{W}$ can exceed 1 even when every $\lvert\rho_{k\ell}\rvert < 1$. The three implementations handle this as follows:

- A evaluates an explicit admissibility statistic, $\min_i\det(\boldsymbol{I}_K - \omega_i\boldsymbol{M})$, computed from the symmetric-polynomial expansion, and treats candidates with a non-positive minimum as non-stationary (graded violation $\geq 1$).

- C rejects any candidate for which some real factor $1 - a_j\omega_i$ falls below $10^{-12}$ (log-likelihood set to $-\infty$); complex factor pairs contribute through their modulus.

- B rejects candidates only where a $2\times2$ factor determinant is exactly singular (log-determinant $-\infty$); it accumulates $\log\lvert\det(\boldsymbol{I}_2 - \omega_i\boldsymbol{R})\rvert$ through LU factorization, i.e. the absolute value, so sign changes of individual factors are not themselves rejected. B's admissible region is therefore slightly larger than C's, and coincides with it on the connected component around $\boldsymbol{R} = \boldsymbol{0}$ in which all real factors stay positive.

Because $\boldsymbol{W}$ in this directory is a row-standardized adjacency matrix derived from a symmetric neighborhood structure, its spectrum is real with $\omega_i \in (-1,1]$ ($\max_i\omega_i = 1$ for row-standardized $\boldsymbol{W}$), and all estimates reported below lie in the interior region where the three admissibility rules agree. For general asymmetric $\boldsymbol{W}$ with complex eigenvalues, the differences above (B: real parts of $\mathrm{eig}(\boldsymbol{W})$ are used; C: complex factors via modulus with a positivity check on real factors; A: real part of the polynomial value) should be reviewed before use. No post-hoc rescaling of $\boldsymbol{R}$ or $\boldsymbol{\Lambda}$ by spectral radius is performed by any implementation.

### 2.4 Precise definition of the effective degrees of freedom

To make Eq. (20) operational, all three implementations use:

- $\boldsymbol{H}$ (denoted $\boldsymbol{I}_{11}$): the numerically evaluated Hessian (finite differences; `numDeriv::hessian` with Richardson extrapolation, or `optimHess` fallback in B) of the negative profile log-likelihood $-\ell_c(\boldsymbol{\theta}_1)$ (step S3 of §4.2) — $\hat{\boldsymbol{\beta}}$ and $\hat{\boldsymbol{\Sigma}}$ are re-profiled by inner GLS at every evaluation point — taken at $\hat{\boldsymbol{\theta}}_1$ in the original $(\boldsymbol{\rho}, \boldsymbol{\lambda})$ parameterization, not in the tanh-transformed space used by B's optimizer. Since effective degrees of freedom are not invariant to nonlinear reparameterization, this choice matters and is stated here explicitly: the reported $\mathrm{deff}$ corresponds to Eq. (20) in the original spatial-parameter space, consistent with Appendix C of the paper. We avoid the term "exact Hessian"; in B's candidate mode, "exact" refers only to evaluating the numerical Hessian rather than bounding it.

- $\boldsymbol{I}$: the identity on the penalized spatial block $\boldsymbol{\theta}_1$ only (dimension = number of free entries of $\boldsymbol{R}$ and $\boldsymbol{\Lambda}$). Regression coefficients $\boldsymbol{\beta}$ (including the temporal AR terms folded into $\boldsymbol{\beta}$ by Eq. (7)) and the free elements of $\boldsymbol{\Sigma}$ are not penalized and contribute their full count:

  $$
  \mathrm{deff}
  = \mathrm{tr}\!\left[\boldsymbol{H}(\boldsymbol{H}+\gamma\boldsymbol{I})^{-1}\right]
  + p + \#\boldsymbol{\Sigma},
  $$

  where $p = \mathrm{ncol}(\text{X-block})$ and $\#\boldsymbol{\Sigma} = K(K+1)/2$ (full $\boldsymbol{\Sigma}$) or $K$ (diagonal $\boldsymbol{\Sigma}$). At $\gamma = 0$, the trace term is set to $\dim(\boldsymbol{\theta}_1)$ exactly.

- Sample size in pBIC:

  $$
  \mathrm{pAIC} = -2\ell + 2\mathrm{deff},
  \qquad
  \mathrm{pBIC} = -2\ell + \log(KnT_{\mathrm{used}})\mathrm{deff},
  $$

  with $K = 2$, $n = 400$, and $T_{\mathrm{used}} = 1$, i.e. $\log(800)$ in both B and C for this dataset. (B: `log(n_eff_total)` with `n_eff_total = 2 * n * TT`; C: `log(dl$K * dl$n)`.)

- $\gamma$-penalty scope: the shrinkage applies only to the six models with a full $\boldsymbol{R}$ or $\boldsymbol{\Lambda}$ block; the five remaining models are evaluated at $\gamma = 0$, so their reported pAIC/pBIC equal ordinary AIC/BIC.

- Negative Hessian eigenvalues: B symmetrizes $\boldsymbol{H}$ and floors negative eigenvalues at zero, issuing a warning whenever an eigenvalue is below $-10^{-6}$; C computes $\mathrm{tr}[\boldsymbol{H}(\boldsymbol{H}+\gamma\boldsymbol{I})^{-1}]$ directly and falls back to the nominal count if the Hessian evaluation fails. Flooring (and the fallback) are numerical approximations, not verified properties: eigenvalues that are materially negative can indicate non-convergence to a local minimum, unstable finite differences, weak identification, or local non-convexity of the profile objective. In the runs reported below, no flooring warning was triggered beyond round-off level, but users should treat any such warning as a diagnostic. Planned diagnostic outputs per fit (see Section 7): minimum eigenvalue, number of negative eigenvalues, largest absolute negative eigenvalue, condition number of $\boldsymbol{H}$, $\mathrm{deff}$ before/after flooring, and a step-size sensitivity check of the finite-difference Hessian.

## 3. Quantitative verification

The implementations were verified quantitatively using the dataset shipped with this directory (`simulated_data_1111_n400_T5.csv`: $n = 400$ spatial units, $K = 2$ responses, $T = 5$ periods; `spatial_weights_n400.csv`: row-standardized $\boldsymbol{W}$). All verification uses the conditional single-period estimation at $t = 2$ described at the top of this document, and — except for Section 3.5 — a single simulated realization; no claims about repeated-sampling performance are made.

**Verification methodology.** Two complementary tools were used. (i) Numerical-equivalence checks (Sections 3.1 and 3.3) were carried out in an independent re-implementation of the four-stage estimation scheme (S1 initialization, inner GLS profile iteration, multi-start quasi-Newton on tanh-transformed variables) with the eigenvalue-wise log-determinant, into which A's symmetric-polynomial log-determinant was substituted as a drop-in replacement; being algebraic in nature, these results are implementation-language independent. (ii) Runtime measurements (Section 3.5) are actual end-to-end wall-clock executions of the three R scripts on a local workstation: R 4.4.2 (`x86_64-w64-mingw32/x64`, `ucrt`) on Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz (4 cores / 8 threads), 64 GB RAM.

### 3.1 Numerical equivalence of the log-likelihood surface

For $K = 2$, the log-determinant schemes of the three implementations are algebraically identical evaluations of the determinant identity of Appendix B:

$$
\det(\boldsymbol{I}_2 - \omega\boldsymbol{M})
= (1 - \omega a_1)(1 - \omega a_2)
= 1 - \omega\mathrm{tr}(\boldsymbol{M}) + \omega^2\det(\boldsymbol{M}),
$$

where $a_1,a_2$ are the eigenvalues of $\boldsymbol{M}$ — the first form is C's, the middle is B's (as a $2\times2$ determinant), and the trace expansion is A's Newton–Girard form (Eqs. (27)–(28)). Numerically, over 5,000 coefficient matrices drawn uniformly from the admissible region, the maximum absolute discrepancy between the symmetric-polynomial and eigenvalue-wise schemes was $5.7\times10^{-14}$ — the level of double-precision round-off — and on the eigenvalues of the actual $\boldsymbol{W}$ both schemes returned $-8.167312779293$ to all displayed digits. On the common admissible region, the three implementations therefore define the same profile log-likelihood surface.

Per-evaluation cost of the two schemes is effectively identical (symmetric-polynomial $\approx 1.1\times$ the eigenvalue-wise scheme in the cross-check implementation), whereas a naive dense $800\times800$ log-determinant that ignores the precomputed eigenvalues of $\boldsymbol{W}$ is roughly three orders of magnitude ($\approx1{,}400\times$) more expensive per call. Since the quasi-Newton loop evaluates the determinant repeatedly, the eigenvalue precomputation of Appendix B of the paper is the dominant efficiency factor at this problem size, and all three implementations adopt it; only the ratios, not the absolute per-call times, are meaningful here.

### 3.2 Full eleven-model estimation ($\gamma = 0$, single realization)

All eleven models were estimated at $\gamma = 0$ through the unified verification pipeline (multi-start optimization included). In this simulation realization, AIC selected the data-generating model (ID 1111, MGNST; 448.06), whereas BIC selected the more parsimonious 1011 (MSAR; 548.78) — a pattern consistent with BIC's heavier complexity penalty, but a single realization neither validates nor refutes the selection performance of either criterion (see Section 7 for the planned Monte Carlo study). Estimates for the 1111 fit were

$$
\hat{\boldsymbol{R}} =
\begin{bmatrix}
0.371 & -0.131\\
0.059 & 0.347
\end{bmatrix},
\qquad
\hat{\boldsymbol{\Lambda}} =
\begin{bmatrix}
0.418 & -0.360\\
0.209 & 0.212
\end{bmatrix},
\qquad
\hat{\boldsymbol{\Sigma}} =
\begin{bmatrix}
0.099 & 0.034\\
0.034 & 0.103
\end{bmatrix},
$$

against data-generating values

$$
\boldsymbol{R} =
\begin{bmatrix}
0.45 & 0.15\\
0.10 & 0.35
\end{bmatrix},
\qquad
\boldsymbol{\Lambda} =
\begin{bmatrix}
0.22 & 0.08\\
0.05 & 0.17
\end{bmatrix},
\qquad
\boldsymbol{\Sigma} =
\begin{bmatrix}
0.10 & 0.03\\
0.03 & 0.10
\end{bmatrix}
$$

(Section 5); all twelve regression coefficients (including the temporal AR terms) were stably estimated. Closed-form models (0011/VARX, 000d) required zero likelihood evaluations, and the mask-decoding pipeline processed every model identically. Sampling variability of a single realization fully accounts for deviations of this magnitude from the true values; no bias claim is made or supported here.

### 3.3 Agreement of estimates across log-determinant schemes

Four representative models were re-estimated in the common verification pipeline with A's symmetric-polynomial log-determinant substituted for the eigenvalue-wise routine, holding everything else (data, initial values, optimizer settings) fixed:

| Model | $\max\lvert\Delta(\boldsymbol{\rho},\boldsymbol{\lambda})\rvert$ | $\max\lvert\Delta\boldsymbol{\beta}\rvert$ | $\max\lvert\Delta\boldsymbol{\Sigma}\rvert$ | $\lvert\Delta\log\mathrm{Lik}\rvert$ |
|---|---:|---:|---:|---:|
| 1111 | $5.8\times10^{-8}$ | $1.1\times10^{-8}$ | $2.4\times10^{-10}$ | 0 (all digits) |
| 1011 | $9.1\times10^{-9}$ | $5.8\times10^{-9}$ | $1.6\times10^{-10}$ | 0 (all digits) |
| 0111 | $2.2\times10^{-7}$ | $6.0\times10^{-9}$ | $2.2\times10^{-9}$ | $5.7\times10^{-14}$ |
| dddd | $7.2\times10^{-9}$ | $1.6\times10^{-9}$ | $1.3\times10^{-10}$ | 0 (all digits) |

The discrepancies are numerical noise at the level of the quasi-Newton convergence tolerance. The supported conclusion is deliberately limited: for the four models examined in the common verification pipeline, replacing the eigenvalue-wise log-determinant routine with A's symmetric-polynomial routine produced numerically indistinguishable estimates. This validates the interchangeability of the two log-determinant schemes — the main mathematical difference among the implementations — but it does not establish that the three actual R scripts produce identical end-to-end output for all eleven models and all $\gamma$; that comparison, run from identical inputs and initial values, remains future work (Section 7). One practical caveat: A's initialization depends on `spatialreg`, so end-to-end execution of A is impossible where that package cannot be installed, in contrast to B and C.

### 3.4 $\gamma$ path, warm starts, and effective degrees of freedom (single realization)

For model 1111 the verification pipeline traced the $\gamma$ path $\{0, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 150\}$ (note: the three R scripts default to the wider grid $\{0\}\cup10^{\{-2,-1.5,\ldots,4\}}$, extending to $10^4$). Because a warm start adds the previous solution to the multi-start set rather than replacing it — the same policy as B's — warm starting increased total cost by $\approx18\%$ ($\approx1.18\times$ likelihood evaluations); in this design it is insurance for path continuity and globality, not an accelerator. An ablation showed that restricting the second and later grid points to a three-member warm-start set yields a $\approx3.3\times$ reduction in likelihood evaluations with objective degradation of at most $1.2\times10^{-12}$ — effectively zero — making it a promising refinement for dense grids.

Under the profile-Hessian effective degrees of freedom, the effective dimension of the spatial parameters shrank from 8 at $\gamma = 0$ to 4.6 at $\gamma = 150$, confirming in this simulated-data example that the shrinkage of Eq. (20) propagates into the penalized criteria. pAIC attained an interior minimum near $\gamma \approx 50$ (446.21). pBIC, in contrast, decreased monotonically over the examined path; its minimum was not attained within the examined grid — the smallest reported pBIC (543.28) occurred at the upper boundary $\gamma = 150$, so no optimal $\gamma$ can be concluded for pBIC from this path. Extending the search (the scripts' default grid reaches $10^4$; further extension to, e.g., $10^5$, or a direct comparison with the $\gamma \to \infty$ limiting model in which the penalized spatial block is fully shrunk) is required before interpreting the pBIC-selected $\gamma$.

### 3.5 End-to-end runtime benchmark (measured in R)

Each script was executed end-to-end on the same workstation (R 4.4.2, Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz, 4C/8T, 64 GB RAM) with the same inputs ($n = 400$, $K = 2$, target period $t = 2$, $\gamma$ grid $\{0\}\cup10^{\{-2,-1.5,\ldots,4\}}$). Wall-clock times were recorded by timing-instrumented wrappers of the scripts:

| Script | Run date | Wall-clock time | Relative to A |
|---|---|---:|---:|
| `implement-A.r` | 2026-06-11 | 1,165.4 s (19.4 min) | $1.0\times$ |
| `implement-B.r` | 2026-06-09/10 | 19,827.4 s (5.51 h) | $17.0\times$ |
| `implement-C.r` | 2026-06-11 | 2,915.0 s (48.6 min) | $2.5\times$ |

**Interpretation — with the appropriate caution.** Section 3.1 established that the likelihood kernels are numerically equivalent per evaluation, so the spread in the table reflects how much work each program performs, not how fast its kernel is. The design differences observable in the code are: B performs, for every penalized model and every $\gamma$ grid point, an exhaustive deterministic multi-start (typically 6–10 quasi-Newton trajectories, with retry escalation), and evaluates numerical profile Hessians whose every finite-difference point triggers an inner GLS re-profiling of $(\boldsymbol{\beta}, \boldsymbol{\Sigma})$ — partially mitigated by candidate screening; C runs a single warm-started L-BFGS-B trajectory per $\gamma$ and one numerical Hessian per $\gamma > 0$; A distributes its experiment modules over PSOCK workers and refines $\gamma$ in a two-stage search, and its master process accumulated only 84.3 s of CPU time (64.8 user + 19.6 system) against 1,165.4 s of wall time, i.e. the computation ran almost entirely in parallel workers. These compounded differences — multi-start count, Hessian evaluation protocol, screening, parallel granularity — plausibly account for the $17\times$ and $2.5\times$ ratios, but each run is a single execution of a distinct end-to-end program, so the table should be read as an operational-cost profile of each program's design, not as an isolated measurement of any single mechanism. Per-stage profiling that decomposes the totals (optimization vs. Hessian vs. I/O) is listed as future work. The robustness trade-off should also be noted: B's multi-start protocol purchases protection against local optima that C's single-trajectory design does not have.

## 4. Recommended usage

- **Routine estimation runs and quick reproduction:** `implement-C.r` — the fastest complete pipeline of the two $K = 2$ programs, with the strictest admissibility check, at the cost of single-start optimization (no multi-start protection against local optima).

- **Publication-grade runs, convergence auditing, and $\gamma$-path analysis:** `implement-B.r` — multi-start with reported diagnostics (`n_starts`, convergence), candidate-screened effective-df evaluation whose selected pAIC/pBIC are exact, $\gamma$-path CSV output, and hardened unattended execution.

- **Extending to $K \geq 3$, or independent cross-checking:** `implement-A.r`, whose symmetric-polynomial log-determinant and explicit admissibility statistic constitute the general-$K$ reference.

## 5. Directory contents and data

| File / directory | Description |
|---|---|
| `implementation_A/implement-A.r`, `implementation_B/implement-B.r`, `implementation_C/implement-C.r` | The three implementations (Section 1) |
| `implementation_A/output_n400/` | A's outputs for the benchmark run (comparison tables, γ searches, model-selection summary) |
| `implementation_B/results_n400/` | B's outputs (per-model estimates, γ paths, run summary, diagnostics) |
| `implementation_C/MSTR_Final_*.xlsx`, `.../MSTR_Final_Long_Format.csv`, `.../MSTR_Final_Execution_Time.txt` | C's outputs (summary workbooks, long-format CSV, execution-time log) |
| `calc_environment/` | Screenshots documenting the R version and benchmark hardware (Section 6) |
| `simulated_data_1111_n400_T5.csv` | Simulated panel: 2,000 rows = 400 regions $\times$ 5 periods; columns `region`, `time`, `y1`, `y2`, `x_common1`, `x_common2`, `x_specific1_1`, `x_specific2_1` |
| `spatial_weights_n400.csv` | $400\times400$ spatial weight matrix, shipped already row-standardized (all row sums = 1); scripts re-validate and normalize with tolerance $10^{-8}$ |

**Data-generating process (model ID 1111, MGNST, Eqs. (4)–(5)).** The simulated data were generated with $K = 2$ and

$$
\boldsymbol{R} =
\begin{bmatrix}
0.45 & 0.15\\
0.10 & 0.35
\end{bmatrix},
\qquad
\boldsymbol{\Lambda} =
\begin{bmatrix}
0.22 & 0.08\\
0.05 & 0.17
\end{bmatrix},
$$

$$
\boldsymbol{A} =
\begin{bmatrix}
0.55 & 0.05\\
0.10 & 0.65
\end{bmatrix},
\qquad
\boldsymbol{\Sigma} =
\begin{bmatrix}
0.10 & 0.03\\
0.03 & 0.10
\end{bmatrix}.
$$

The regression coefficients are

$$
\boldsymbol{\beta}:
\quad
\begin{aligned}
\text{intercepts} &: (0.20,-0.10),\\
x_{\mathrm{common1}} &: (1.20,0.80),\\
x_{\mathrm{common2}} &: (-0.60,-0.30),\\
x_{\mathrm{specific1\_1}} &: 0.80,\\
x_{\mathrm{specific2\_1}} &: 0.50,
\end{aligned}
$$

where the first entry corresponds to equation $y_1$ and the second to $y_2$. These values are recorded as `TRUE_PARAMS` in `implement-A.r`.

**Parameter matrices and dimensions (Eqs. (3)–(4)).** $\boldsymbol{R}$, $\boldsymbol{\Lambda}$, and $\boldsymbol{A}$ are $K\times K = 2\times2$ coefficient matrices for the spatial lag of the responses, the spatial lag of the errors, and the temporal AR(1) terms, respectively; $\boldsymbol{\Sigma}$ is the $2\times2$ error covariance; $\boldsymbol{\beta}$ stacks, per equation, an intercept, two common regressors, one equation-specific regressor, and the AR(1) terms folded in by Eq. (7) — 12 regression-block coefficients in total for the full models.

**Model IDs (Figure 1 of the paper).** Each four-character ID codes $(\boldsymbol{R}, \boldsymbol{\Lambda}, \boldsymbol{A}, \boldsymbol{\Sigma})$ as full $K\times K$ (`1`), diagonal (`d`), or zero (`0`):

| ID | Constraints | Free parameters ($K = 2$, this design) | Name |
|---|---|---:|---|
| 0011 | $\boldsymbol{R} = \boldsymbol{\Lambda} = \boldsymbol{O}$ | $\boldsymbol{\beta}(12) + \boldsymbol{\Sigma}(3) = 15$ | VARX |
| 1011 | $\boldsymbol{\Lambda} = \boldsymbol{O}$ | $+\boldsymbol{R}(4) = 19$ | MSAR |
| 0111 | $\boldsymbol{R} = \boldsymbol{O}$ | $+\boldsymbol{\Lambda}(4) = 19$ | MSEM |
| 1111 | — | $+\boldsymbol{R}(4) + \boldsymbol{\Lambda}(4) = 23$ | MGNST |
| 1001 | $\boldsymbol{\Lambda} = \boldsymbol{A} = \boldsymbol{O}$ | $\boldsymbol{\beta}(8) + \boldsymbol{R}(4) + \boldsymbol{\Sigma}(3) = 15$ | MSAR w/o time-AR |
| 0101 | $\boldsymbol{R} = \boldsymbol{A} = \boldsymbol{O}$ | $\boldsymbol{\beta}(8) + \boldsymbol{\Lambda}(4) + \boldsymbol{\Sigma}(3) = 15$ | MSEM w/o time-AR |
| 1101 | $\boldsymbol{A} = \boldsymbol{O}$ | $\boldsymbol{\beta}(8) + \boldsymbol{R}(4) + \boldsymbol{\Lambda}(4) + \boldsymbol{\Sigma}(3) = 19$ | MGNS w/o time-AR |
| 000d | $\boldsymbol{R} = \boldsymbol{\Lambda} = \boldsymbol{A} = \boldsymbol{O}$, $\boldsymbol{\Sigma}$ diagonal | $\boldsymbol{\beta}(8) + \boldsymbol{\Sigma}(2) = 10$ | Independent regressions |
| d0dd | $\boldsymbol{R}$, $\boldsymbol{A}$, and $\boldsymbol{\Sigma}$ diagonal, $\boldsymbol{\Lambda} = \boldsymbol{O}$ | $\boldsymbol{\beta}(8+2) + \boldsymbol{R}(2) + \boldsymbol{\Sigma}(2) = 14$ | Independent SAR |
| 0ddd | $\boldsymbol{\Lambda}$, $\boldsymbol{A}$, and $\boldsymbol{\Sigma}$ diagonal, $\boldsymbol{R} = \boldsymbol{O}$ | $\boldsymbol{\beta}(8+2) + \boldsymbol{\Lambda}(2) + \boldsymbol{\Sigma}(2) = 14$ | Independent SEM |
| dddd | all diagonal | $\boldsymbol{\beta}(8+2) + \boldsymbol{R}(2) + \boldsymbol{\Lambda}(2) + \boldsymbol{\Sigma}(2) = 16$ | Independent GNS |

$\boldsymbol{\beta}$ counts include intercepts and, where $\boldsymbol{A} \neq \boldsymbol{O}$, the AR terms: 2 per equation for `1`, 1 per equation for `d`.

## 6. Reproducibility

**Environment used for all reported runs.** R 4.4.2 (2024-10-31 ucrt, "Pile of Leaves"), platform `x86_64-w64-mingw32/x64`, Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz (4 physical / 8 logical cores), 64 GB RAM, reference BLAS as shipped with the Windows R binary. The R version and hardware are documented by the screenshots in `calc_environment/`; users are encouraged to pin package versions with `renv::snapshot()`.

**Determinism.** No random number generation is used in estimation: A initializes from `spatialreg` fits, B constructs a deterministic multi-start set, and C starts from the constant vector 0.01 — so no seed is required and reruns are bit-reproducible up to BLAS/thread scheduling effects. Optimizer settings: B — `optim(method = "BFGS")` in tanh space, `reltol = 1e-8`, `maxit = 220` (retry 320); C — `optim(method = "L-BFGS-B")`, bounds $\pm0.99$, `maxit = 500`; inner GLS: 50 iterations, tolerance $10^{-6}$ (both). $\boldsymbol{\Sigma}$ is symmetrized as $(\boldsymbol{S} + \boldsymbol{S}^{\mathsf{T}})/2$ at every update; C additionally floors $\boldsymbol{\Sigma}$'s eigenvalues at $10^{-8}$ (`regularize_symmetric`), and B rejects any GLS iterate whose $\boldsymbol{\Sigma}$ has a non-positive eigenvalue.

**Parallel configuration on the benchmark machine.** A: PSOCK cluster, default `MSTR_CORES` = physical cores $- 1 = 3$ (overridable via env var `MSTR_CORES`), BLAS threads pinned to 1 per worker; B: per-model PSOCK, auto rule for $n = 400$ gives $\lfloor0.75(8 - 1)\rfloor = 5$ workers; C: `doParallel`, $\min(11, 8 - 1) = 7$ workers.

**Input file integrity (SHA-256).**

```text
7c0f63b34ab6d491241b7e277808e289f09f3d188026b0a94b523bc6eafc1fac  simulated_data_1111_n400_T5.csv
46269338737e6cbcd37ed60084b72793da6bbaee7538987b394f3390be4d235d  spatial_weights_n400.csv
98bb5b2c355849cecd5f9137adab25037700060da60195926aa1309d7dac464b  implementation_A/implement-A.r
764fdf056d89949b11099c7a9fa3248d99e9535881f3f48383e487fd80e48cfe  implementation_B/implement-B.r
dbbdfd44aa8e90d86170aee88b399558dcce74f792639a72b185d1c379a89ec7  implementation_C/implement-C.r
```

**How to run.** All three scripts resolve the data files relative to the script directory as `../simulated_data_1111_n400_T5.csv` and `../spatial_weights_n400.csv`; place each script one directory below the data (or edit `DATA_FILE`/`data_file` at the top of the script).

```bash
# implement-B.r: non-interactive runs must set the output root policy
MSTR_USE_CURRENT_DIR=yes Rscript implement-B.r  # or MSTR_OUT_ROOT=/path/to/out

# implement-C.r
Rscript implement-C.r

# implement-A.r (requires spdep/spatialreg; worker count via MSTR_CORES)
MSTR_CORES=4 Rscript implement-A.r
```

**Criteria formulas as computed (for this dataset):**

$$
\mathrm{pAIC} = -2\ell(\hat{\boldsymbol{\theta}}) + 2\mathrm{deff},
\qquad
\mathrm{pBIC} = -2\ell(\hat{\boldsymbol{\theta}}) + \log(800)\mathrm{deff},
$$

with $\mathrm{deff}$ as defined in Section 2.4; ordinary AIC/BIC use the nominal parameter counts in the table in Section 5. The averaged pseudo-$\bar{R}^2$ (§4.4 of the paper) is reported descriptively by B.

## 7. Limitations and future work

1. The runtime benchmark of Section 3.5 consists of a single end-to-end run of each script on one Windows workstation; replicated timings and per-stage profiling (optimization vs. Hessian evaluation vs. I/O) are needed before attributing the observed ratios to specific mechanisms.

2. All statistical verification rests on a single simulated realization and on the conditional single-period estimation at $t = 2$; a Monte Carlo study ($\geq 500$ replications at $n = 100/400/900$) reporting true-model selection rates, over-/under-selection rates, parameter bias and RMSE, convergence rates, and the distribution of the selected $\gamma$ is required before any claim about the selection performance of pAIC/pBIC, as is an extension to joint estimation over all $T$ periods.

3. The cross-implementation agreement of Section 3.3 covers four models within a common pipeline; a full end-to-end comparison of the three actual R scripts — all eleven models, the entire $\gamma$ grid, identical inputs and initial values — remains to be run and published.

4. The pBIC path of Section 3.4 terminated at the boundary of the examined grid; the grid must be extended (or the $\gamma \to \infty$ limit compared directly) before the pBIC-selected $\gamma$ is interpreted.

5. The Hessian-eigenvalue flooring and fallback rules of Section 2.4 are approximations; the per-fit diagnostics listed there (minimum eigenvalue, negative-eigenvalue count, condition number, $\mathrm{deff}$ before/after flooring, finite-difference step sensitivity) should be added to the standard output.

6. Porting A's symmetric-polynomial log-determinant into B/C to lift the $K = 2$ restriction, and unifying the three programs' admissibility rules (Section 2.3), are natural next steps.
