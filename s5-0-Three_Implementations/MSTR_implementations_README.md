# Three Reference Implementations of Simultaneous Parameter Estimation for Multivariate Spatio-Temporal Regression (MSTR)

This repository provides three **separately developed, related R implementations** — `implement-A.r`, `implement-B.r`, and `implement-C.r` — of the maximum likelihood estimation procedure for the multivariate general nesting spatial (MGNS) model and its ten nested submodels, as formulated in the accompanying paper (`estimation-algorithm-mstr.pdf`; Eqs. (8), (13)–(16), (17)–(23), (24)–(26), and Appendix A). All three follow the four-stage scheme of Subsection 3.3 of the paper: (S1) univariate-GNS initial values, (S2) inner iterative GLS updates of β̂ (Eq. 13) and Σ̂ (Eq. 15), (S3) profile likelihood in the spatial parameters (Eq. 35), and (S4) quasi-Newton updates of the free entries of (R, Λ), with penalized information criteria pAIC/pBIC (Eq. 25) based on the effective degrees of freedom of Eq. (24).

The three programs are not fully independent — they implement the same specification and share design conventions — but they were written separately and differ in code base, log-determinant scheme, optimization strategy, and operational design. **A provides a partially independent numerical reference** (general-K symmetric-function log-determinant, external `spatialreg` initialization), while **B and C are two architecturally distinct realizations of the K = 2 estimator**. Agreement among their outputs therefore constitutes a meaningful, though not fully independent, cross-check.

> **Scope of the numerical verification.** All verification and benchmark results below concern the **conditional estimation at a single target period t = 2**, with the t = 1 observations entering only through the AR(1) regressors of Eq. (7). They do **not** verify a joint likelihood over all five periods of the supplied panel. (The dataset contains T = 5 periods; `TARGET_TIME`/`time_now` is set to 2 in all three scripts.)

## 1. The three implementations

### `implement-A.r` — Reference implementation (general K)

A comprehensive development version (13,994 lines) integrating all modules produced during development: data preparation, eleven model-specific fitters, penalized spatial estimation, CSV output, and a parallel experiment runner. Its distinguishing feature is a faithful implementation of the log-determinant log|I − R⊗W| via the elementary-symmetric-polynomial expansion of Eqs. (18)–(23), computed from power-sum traces r_j = tr(Rʲ) by Newton's identities and valid for general response dimension K ≤ 7. Admissibility of candidate (R, Λ) is checked explicitly by a `spectral_radius()` routine that tests min_i det(I_K − ω_i M) > 0 over all eigenvalues ω_i of W and reports a graded violation measure. Initial values are obtained through `spatialreg::lagsarlm` / `errorsarlm`, so A carries the largest dependency footprint (`spdep`, `spatialreg`, `Matrix`, `numDeriv`, `parallel`). The γ search is two-stage (coarse grid, then log-scale refinement around the coarse optimum). A is retained as (i) the reference for extending the estimator beyond K = 2 and (ii) a partially independent implementation for cross-checking B and C.

### `implement-B.r` — Fully instrumented unified runner (K = 2)

A self-contained runner (2,330 lines) specialized to K = 2. All eleven model IDs are decoded into masks on (R, Λ, A, Σ) and estimated through a single shared pipeline. Optimization is by multi-start BFGS (`optim`, `reltol = 1e-8`, `maxit` 220 with a 320-iteration retry) on **elementwise tanh-transformed** variables, with a deterministic multi-start set (univariate initial values, zero starts, sign and scale perturbations) augmented — not replaced — by warm starts along the γ path. The log-determinant is evaluated eigenvalue-wise as Σ_i log|det(I₂ − ω_i R)|, the K = 2 form of Eq. (18). Effective degrees of freedom use a numerically evaluated profile Hessian (Section 2.4 below) with a **"candidate" screening mode** (`DEFF_EVALUATION_MODE = "candidate"`, the default): rigorous lower/upper pAIC/pBIC bounds identify the γ grid rows that could still be optimal, and the expensive numerical Hessian is evaluated only for those rows, so the *selected* pAIC/pBIC values are the same as under exhaustive evaluation. Every fit reports `n_starts` and `convergence` diagnostics. B also contains the operational hardening previously attributed elsewhere: interactive/environment-variable control of the output root (`MSTR_OUT_ROOT`, `MSTR_USE_CURRENT_DIR`), and safe CSV/workbook writes that fall back to timestamped filenames when files are locked. Dependencies: `Matrix`, `openxlsx`, `parallel` (`numDeriv` optional; `optimHess` fallback).

### `implement-C.r` — Compact direct implementation (K = 2)

C (646 lines) is a compact, direct transcription of the estimation scheme. The profile likelihood forms I_Kn − R⊗W and I_Kn − Λ⊗W explicitly (dense 800×800 for this dataset) for the GLS steps, while the log-determinant is evaluated from the eigenvalues of the K × K coefficient matrix and of W as Σ_{i,j} log|1 − a_j ω_i| with a strict admissibility check (any real factor ≤ 10⁻¹² ⇒ the candidate is rejected with objective 10¹⁰). Optimization is single-trajectory `optim(method = "L-BFGS-B")` with box constraints −0.99 ≤ θ ≤ 0.99 on each free spatial parameter and warm starts carried along the γ grid (initial value 0.01 for every free parameter at the first grid point; no multi-start). Effective degrees of freedom use a `numDeriv::hessian` of the negative profile log-likelihood at each γ > 0, with a documented fallback to the nominal spatial parameter count when the Hessian computation fails. Models are parallelized with `doParallel`/`foreach` (one model per worker), and results are written to xlsx. Dependencies: `Matrix`, `parallel`, `doParallel`, `foreach`, `openxlsx`, `numDeriv`.

> **Correction relative to an earlier draft of this document.** An earlier draft stated that B and C share a byte-identical numerical core differing only in I/O, and attributed candidate-γ screening to C. Inspection of the published files shows the opposite: **candidate screening is implemented in B and absent from C**, and B and C are architecturally distinct programs (multi-start tanh-space BFGS with screening in B; single-start box-constrained L-BFGS-B with dense Kronecker algebra in C). All statements below are grounded in the published source files; their SHA-256 checksums are listed in Section 6.

## 2. Design comparison

| Criterion | A | B | C |
| --- | --- | --- | --- |
| Size / structure | 13,994 lines, all development modules | 2,330 lines, unified runner | 646 lines, compact direct implementation |
| log-det of Eq. (18) | Symmetric-polynomial expansion via Newton's identities, general K ≤ 7 | Eigenvalue-wise 2×2 determinants Σ log\|det(I₂ − ω_i R)\| | Eigenvalue-pair products Σ log\|1 − a_j ω_i\| with strict positivity check |
| Constraint handling (S4) | Explicit admissibility check min_i det(I_K − ω_i M) > 0 | Elementwise tanh (\|ρ_kℓ\|, \|λ_kℓ\| < 1); singular determinants rejected (−∞) | Box constraints ±0.99; real factors ≤ 10⁻¹² rejected (−∞) |
| Optimization | `optim` per module (incl. penalized L-BFGS-B); numerical Hessians | Multi-start BFGS in tanh space, warm starts add to start set, retry logic | Single warm-started L-BFGS-B trajectory per γ |
| Effective d.o.f. (Eq. 24) | tr[I₁₁(I₁₁+γI)⁻¹] + p + K(K+1)/2, numerical Hessian (GIC module) | Numerical profile θ-Hessian; candidate-γ screening with exact-at-selection guarantee; eigenvalues floored at 0 with warning | Numerical profile θ-Hessian per γ > 0; fallback to nominal count on failure |
| γ grid | {0} ∪ 10^{−2:4 by 0.5} coarse, then log-scale refinement | {0} ∪ 10^{−2:4 by 0.5} (14 points) | {0} ∪ 10^{−2:4 by 0.5} (14 points) |
| γ-penalized models | Full spatial-block models | The six IDs {1111, 0111, 1011, 1101, 0101, 1001}; others at γ = 0 with ordinary AIC/BIC | Same six (R- or Λ-code = "1"); others at γ = 0 |
| Initial values (S1) | `spatialreg::lagsarlm` / `errorsarlm` | In-house univariate GNS (cached) | Constant 0.01 per free spatial parameter |
| Parallelism | PSOCK cluster, `MSTR_CORES` (default: physical cores − 1), BLAS threads pinned to 1 | Per-model PSOCK; auto worker count (for n = 400: ⌊0.75·(logical cores − 1)⌋), sequential retry | `doParallel`, min(11, logical cores − 1) workers |
| Output | Model-specific CSV writers | xlsx + run-summary / γ-path CSV, safe writes, output-root prompt/env vars | xlsx summary tables |
| Dependencies | many (`spdep`, `spatialreg`, `systemfit`, …) | few (`Matrix`, `openxlsx`, `parallel`; `numDeriv` optional) | moderate (`Matrix`, `doParallel`, `foreach`, `openxlsx`, `numDeriv`) |

### 2.3 What is actually constrained (stationarity and admissibility)

The paper's Step S4 states the elementwise constraints |ρ_kℓ| < 1, |λ_kℓ| < 1. As implemented, these elementwise bounds are enforced by tanh transformation (B) or L-BFGS-B boxes (C, A's penalized modules). **Elementwise bounds alone do not guarantee stationarity/invertibility of the matrix process**: for non-diagonal R, the relevant condition is that all factors 1 − a_j(R)·ω_i remain bounded away from zero (equivalently det(I_Kn − R⊗W) ≠ 0 along a path connected to R = 0), and the spectral radius of R⊗W can exceed 1 even when every |ρ_kℓ| < 1. The three implementations handle this as follows:

- **A** evaluates an explicit admissibility statistic, min_i det(I_K − ω_i M) computed from the symmetric-polynomial expansion, and treats candidates with a non-positive minimum as non-stationary (graded violation ≥ 1).
- **C** rejects any candidate for which some real factor 1 − a_j ω_i falls below 10⁻¹² (log-likelihood set to −∞); complex factor pairs contribute through their modulus.
- **B** rejects candidates only where a 2×2 factor determinant is exactly singular (log-determinant −∞); it accumulates log|det(I₂ − ω_i R)| through LU factorization, i.e. the absolute value, so sign changes of individual factors are not themselves rejected. B's admissible region is therefore slightly larger than C's, and coincides with it on the connected component around R = 0 in which all real factors stay positive.

Because W in this repository is a row-standardized adjacency matrix derived from a symmetric neighborhood structure, its spectrum is real with ω_i ∈ (−1, 1] (max ω_i = 1 for row-standardized W), and all estimates reported below lie in the interior region where the three admissibility rules agree. For general asymmetric W with complex eigenvalues, the differences above (B: real parts of eig(W) are used; C: complex factors via modulus with a positivity check on real factors; A: real part of the polynomial value) should be reviewed before use. No post-hoc rescaling of R or Λ by spectral radius is performed by any implementation.

### 2.4 Precise definition of the effective degrees of freedom

To make Eq. (24) operational, all three implementations use:

- **H (denoted I₁₁):** the **numerically evaluated** Hessian (finite differences; `numDeriv::hessian` with Richardson extrapolation, or `optimHess` fallback in B) of the **negative profile log-likelihood** −ℓ_c(θ₁) of Eq. (35) — β̂ and Σ̂ are re-profiled by inner GLS at every evaluation point — taken **at θ̂₁ in the original (ρ, λ) parameterization**, *not* in the tanh-transformed space used by B's optimizer. Since effective degrees of freedom are not invariant to nonlinear reparameterization, this choice matters and is stated here explicitly: the reported deff corresponds to Eq. (24) in the original spatial-parameter space, consistent with Appendix A of the paper. We avoid the term "exact Hessian"; in B's candidate mode, "exact" refers only to *evaluating* the numerical Hessian rather than bounding it.
- **I:** the identity on the penalized spatial block θ₁ only (dimension = number of free entries of R and Λ). Regression coefficients β (including the temporal AR terms folded into β by Eq. (7)) and the free elements of Σ are **not** penalized and contribute their full count: deff = tr[H(H + γI)⁻¹] + p + #Σ, where p = ncol(X-block) and #Σ = K(K+1)/2 (full Σ) or K (diagonal Σ). At γ = 0 the trace term is set to dim(θ₁) exactly (Eq. 44).
- **Sample size in pBIC:** pAIC = −2ℓ + 2·deff and pBIC = −2ℓ + log(Kn·T_used)·deff with K = 2, n = 400, T_used = 1, i.e. **log(800)** in both B and C for this dataset. (B: `log(n_eff_total)` with `n_eff_total = 2·n·TT`; C: `log(dl$K * dl$n)`.)
- **γ-penalty scope:** the shrinkage applies only to the six models with a full R or Λ block; the five remaining models are evaluated at γ = 0, so their reported pAIC/pBIC equal ordinary AIC/BIC.

**Negative Hessian eigenvalues.** B symmetrizes H and floors negative eigenvalues at zero, issuing a warning whenever an eigenvalue is below −10⁻⁶; C computes tr[H(H+γI)⁻¹] directly and falls back to the nominal count if the Hessian evaluation fails. Flooring (and the fallback) are **numerical approximations, not verified properties**: eigenvalues that are materially negative can indicate non-convergence to a local minimum, unstable finite differences, weak identification, or local non-convexity of the profile objective. In the runs reported below, no flooring warning was triggered beyond round-off level, but users should treat any such warning as a diagnostic. Planned diagnostic outputs per fit (see Section 7): minimum eigenvalue, number of negative eigenvalues, largest |negative| eigenvalue, condition number of H, deff before/after flooring, and a step-size sensitivity check of the finite-difference Hessian.

## 3. Quantitative verification

The implementations were verified quantitatively using the dataset shipped with this repository (`simulated_data_1111_n400_T5.csv`: n = 400 spatial units, K = 2 responses, T = 5 periods; `spatial_weights_n400.csv`: row-standardized W). **All verification uses the conditional single-period estimation at t = 2** described at the top of this document, and — except for Section 3.5 — a **single simulated realization**; no claims about repeated-sampling performance are made.

**Verification methodology.** Two complementary tools were used. (i) *Numerical-equivalence checks* (Sections 3.1 and 3.3) were carried out in an independent re-implementation of the four-stage estimation scheme (S1 initialization, inner GLS profile iteration, multi-start quasi-Newton on tanh-transformed variables) with the eigenvalue-wise log-determinant, into which A's symmetric-polynomial log-determinant was substituted as a drop-in replacement; being algebraic in nature, these results are implementation-language independent. (ii) *Runtime measurements* (Section 3.5) are actual end-to-end wall-clock executions of the three R scripts on a local workstation: R 4.4.2 (x86_64-w64-mingw32/x64, ucrt) on Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz (4 cores / 8 threads), 64 GB RAM.

### 3.1 Numerical equivalence of the log-likelihood surface

For K = 2 the log-determinant schemes of the three implementations are algebraically identical evaluations of Eq. (18): det(I₂ − ωM) = (1 − ω a₁)(1 − ω a₂) = 1 − ω·tr M + ω²·det M, where a₁, a₂ are the eigenvalues of M — the first form is C's, the middle is B's (as a 2×2 determinant), and the trace expansion is A's Eqs. (19)–(20). Numerically, over 5,000 coefficient matrices drawn uniformly from the admissible region, the maximum absolute discrepancy between the symmetric-polynomial and eigenvalue-wise schemes was 5.7×10⁻¹⁴ — the level of double-precision round-off — and on the eigenvalues of the actual W both schemes returned −8.167312779293 to all displayed digits. On the common admissible region, the three implementations therefore define the same profile log-likelihood surface.

Per-evaluation cost of the two schemes is effectively identical (symmetric-polynomial ≈ 1.1× the eigenvalue-wise scheme in the cross-check implementation), whereas a naive dense 800×800 log-determinant that ignores the precomputed eigenvalues of W is roughly three orders of magnitude (≈1,400×) more expensive per call. Since the quasi-Newton loop evaluates the determinant repeatedly, the eigenvalue precomputation of Subsection 3.4 of the paper is the dominant efficiency factor at this problem size, and all three implementations adopt it; only the ratios, not the absolute per-call times, are meaningful here.

### 3.2 Full eleven-model estimation (γ = 0, single realization)

All eleven models were estimated at γ = 0 through the unified verification pipeline (multi-start optimization included). **In this simulation realization**, AIC selected the data-generating model (ID 1111, MGNS; 448.06), whereas BIC selected the more parsimonious 1011 (MSAR; 548.78) — a pattern consistent with BIC's heavier complexity penalty, but a single realization neither validates nor refutes the selection performance of either criterion (see Section 7 for the planned Monte Carlo study). Estimates for the 1111 fit were R̂ = [0.371, −0.131; 0.059, 0.347], Λ̂ = [0.418, −0.360; 0.209, 0.212], Σ̂ = [0.099, 0.034; 0.034, 0.103], against data-generating values R = [0.45, 0.15; 0.10, 0.35], Λ = [0.22, 0.08; 0.05, 0.17], Σ = [0.10, 0.03; 0.03, 0.10] (Section 5); all twelve regression coefficients (including the temporal AR terms) were stably estimated. Closed-form models (0011/VARX, 000d) required zero likelihood evaluations, and the mask-decoding pipeline processed every model identically. Sampling variability of a single realization fully accounts for deviations of this magnitude from the true values; no bias claim is made or supported here.

### 3.3 Agreement of estimates across log-determinant schemes

Four representative models were re-estimated in the common verification pipeline with A's symmetric-polynomial log-determinant substituted for the eigenvalue-wise routine, holding everything else (data, initial values, optimizer settings) fixed:

| Model | max\|Δ(ρ, λ)\| | max\|Δβ\| | max\|ΔΣ\| | \|Δ logLik\| |
| --- | --- | --- | --- | --- |
| 1111 | 5.8×10⁻⁸ | 1.1×10⁻⁸ | 2.4×10⁻¹⁰ | 0 (all digits) |
| 1011 | 9.1×10⁻⁹ | 5.8×10⁻⁹ | 1.6×10⁻¹⁰ | 0 (all digits) |
| 0111 | 2.2×10⁻⁷ | 6.0×10⁻⁹ | 2.2×10⁻⁹ | 5.7×10⁻¹⁴ |
| dddd | 7.2×10⁻⁹ | 1.6×10⁻⁹ | 1.3×10⁻¹⁰ | 0 (all digits) |

The discrepancies are numerical noise at the level of the quasi-Newton convergence tolerance. The supported conclusion is deliberately limited: **for the four models examined in the common verification pipeline, replacing the eigenvalue-wise log-determinant routine with A's symmetric-polynomial routine produced numerically indistinguishable estimates.** This validates the interchangeability of the two log-determinant schemes — the main mathematical difference among the implementations — but it does **not** establish that the three actual R scripts produce identical end-to-end output for all eleven models and all γ; that comparison, run from identical inputs and initial values, remains future work (Section 7). One practical caveat: A's initialization depends on `spatialreg`, so end-to-end execution of A is impossible where that package cannot be installed, in contrast to B and C.

### 3.4 γ path, warm starts, and effective degrees of freedom (single realization)

For model 1111 the verification pipeline traced the γ path {0, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 150} (note: the three R scripts default to the wider grid {0} ∪ 10^{−2:4 by 0.5}, extending to 10⁴). Because a warm start *adds* the previous solution to the multi-start set rather than replacing it — the same policy as B's — warm starting increased total cost by ≈18% (≈1.18× likelihood evaluations); in this design it is insurance for path continuity and globality, not an accelerator. An ablation showed that restricting the second and later grid points to a three-member warm-start set yields a ≈3.3× reduction in likelihood evaluations with objective degradation of at most 1.2×10⁻¹² — effectively zero — making it a promising refinement for dense grids.

Under the profile-Hessian effective degrees of freedom, the effective dimension of the spatial parameters shrank from 8 at γ = 0 to 4.6 at γ = 150, confirming **in this simulated-data example** that the shrinkage of Eq. (24) propagates into the penalized criteria. pAIC attained an interior minimum near γ ≈ 50 (446.21). pBIC, in contrast, decreased monotonically over the examined path; **its minimum was not attained within the examined grid — the smallest reported pBIC (543.28) occurred at the upper boundary γ = 150**, so no optimal γ can be concluded for pBIC from this path. Extending the search (the scripts' default grid reaches 10⁴; further extension to, e.g., 10⁵, or a direct comparison with the γ → ∞ limiting model in which the penalized spatial block is fully shrunk) is required before interpreting the pBIC-selected γ.

### 3.5 End-to-end runtime benchmark (measured in R)

Each script was executed end-to-end on the same workstation (R 4.4.2, Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz, 4C/8T, 64 GB RAM) with the same inputs (n = 400, K = 2, target period t = 2, γ grid {0} ∪ 10^{−2:4 by 0.5}). Wall-clock times were recorded by timing-instrumented wrappers of the scripts:

| Script | Run date | Wall-clock time | Relative to A |
| --- | --- | --- | --- |
| `implement-A.r` | 2026-06-11 | 1,165.4 s (19.4 min) | 1.0× |
| `implement-B.r` | 2026-06-09/10 | 19,827.4 s (5.51 h) | 17.0× |
| `implement-C.r` | 2026-06-11 | 2,915.0 s (48.6 min) | 2.5× |

**Interpretation — with the appropriate caution.** Section 3.1 established that the likelihood kernels are numerically equivalent per evaluation, so the spread in the table reflects how much work each program *performs*, not how fast its kernel is. The design differences observable in the code are: B performs, for every penalized model and every γ grid point, an exhaustive deterministic multi-start (typically 6–10 quasi-Newton trajectories, with retry escalation), and evaluates numerical profile Hessians whose every finite-difference point triggers an inner GLS re-profiling of (β, Σ) — partially mitigated by candidate screening; C runs a single warm-started L-BFGS-B trajectory per γ and one numerical Hessian per γ > 0; A distributes its experiment modules over PSOCK workers and refines γ in a two-stage search, and its master process accumulated only 84.3 s of CPU time (64.8 user + 19.6 system) against 1,165.4 s of wall time, i.e. the computation ran almost entirely in parallel workers. These compounded differences — multi-start count, Hessian evaluation protocol, screening, parallel granularity — plausibly account for the 17× and 2.5× ratios, but **each run is a single execution of a distinct end-to-end program, so the table should be read as an operational-cost profile of each program's design, not as an isolated measurement of any single mechanism**. Per-stage profiling that decomposes the totals (optimization vs. Hessian vs. I/O) is listed as future work. The robustness trade-off should also be noted: B's multi-start protocol purchases protection against local optima that C's single-trajectory design does not have.

## 4. Recommended usage

- **Routine estimation runs and quick reproduction:** `implement-C.r` — the fastest complete pipeline of the two K = 2 programs, with the strictest admissibility check, at the cost of single-start optimization (no multi-start protection against local optima).
- **Publication-grade runs, convergence auditing, and γ-path analysis:** `implement-B.r` — multi-start with reported diagnostics (`n_starts`, `convergence`), candidate-screened effective-df evaluation whose selected pAIC/pBIC are exact, γ-path CSV output, and hardened unattended execution.
- **Extending to K ≥ 3, or independent cross-checking:** `implement-A.r`, whose symmetric-polynomial log-determinant and explicit admissibility statistic constitute the general-K reference.

## 5. Repository contents and data

| File | Description |
| --- | --- |
| `estimation-algorithm-mstr.pdf` | The accompanying paper (model, estimation scheme, Eq. numbers cited here) |
| `implement-A.r`, `implement-B.r`, `implement-C.r` | The three implementations (Section 1) |
| `simulated_data_1111_n400_T5.csv` | Simulated panel: 2,000 rows = 400 regions × 5 periods; columns `region`, `time`, `y1`, `y2`, `x_common1`, `x_common2`, `x_specific1_1`, `x_specific2_1` |
| `spatial_weights_n400.csv` | 400 × 400 spatial weight matrix, shipped **already row-standardized** (all row sums = 1); scripts re-validate and normalize with tolerance 10⁻⁸ |

**Data-generating process (model ID 1111, MGNS, Eq. (8)).** The simulated data were generated with K = 2 and
R = [0.45, 0.15; 0.10, 0.35], Λ = [0.22, 0.08; 0.05, 0.17], A = [0.55, 0.05; 0.10, 0.65], Σ = [0.10, 0.03; 0.03, 0.10],
β: intercepts (0.20, −0.10), `x_common1` (1.20, 0.80), `x_common2` (−0.60, −0.30), `x_specific1_1` = 0.80, `x_specific2_1` = 0.50 (first entry: equation y₁, second: y₂). These values are recorded as `TRUE_PARAMS` in `implement-A.r`.

**Parameter matrices and dimensions (Eq. (4)).** R, Λ, A are K × K = 2 × 2 coefficient matrices for the spatial lag of the responses, the spatial lag of the errors, and the temporal AR(1) terms respectively; Σ is the 2 × 2 error covariance; β stacks, per equation, an intercept, two common regressors, one equation-specific regressor, and the AR(1) terms folded in by Eq. (7) — 12 regression-block coefficients in total for the full models.

**Model IDs (Table 1 of the paper).** Each four-character ID codes (R, Λ, A, Σ) as full K × K (`1`), diagonal (`d`), or zero (`0`):

| ID | Constraints | Free parameters (K = 2, this design) | Name |
| --- | --- | --- | --- |
| 0011 | R = Λ = O | β(12) + Σ(3) = 15 | VARX |
| 1011 | Λ = O | + R(4) = 19 | MSAR |
| 0111 | R = O | + Λ(4) = 19 | MSEM |
| 1111 | — | + R(4) + Λ(4) = 23 | MGNS |
| 1001 | Λ = A = O | β(8) + R(4) + Σ(3) = 15 | MSAR w/o time-AR |
| 0101 | R = A = O | β(8) + Λ(4) + Σ(3) = 15 | MSEM w/o time-AR |
| 1101 | A = O | β(8) + R(4) + Λ(4) + Σ(3) = 19 | MGNS w/o time-AR |
| 000d | R = Λ = A = O, Σ diag | β(8) + Σ(2) = 10 | Independent regressions |
| d0dd | R, A, Σ diag, Λ = O | β(8+2) + R(2) + Σ(2) = 14 | Independent SAR |
| 0ddd | Λ, A, Σ diag, R = O | β(8+2) + Λ(2) + Σ(2) = 14 | Independent SEM |
| dddd | all diagonal | β(8+2) + R(2) + Λ(2) + Σ(2) = 16 | Independent GNS |

(β counts include intercepts and, where A ≠ O, the AR terms: 2 per equation for `1`, 1 per equation for `d`.)

## 6. Reproducibility

**Environment used for all reported runs.** R 4.4.2 (2024-10-31 ucrt, "Pile of Leaves"), platform x86_64-w64-mingw32/x64, Windows 11 Home, Intel Core i7-6700K @ 4.00 GHz (4 physical / 8 logical cores), 64 GB RAM, reference BLAS as shipped with the Windows R binary. `sessionInfo()` output and package versions for the benchmark runs will be committed alongside a release tag; users are encouraged to pin versions with `renv::snapshot()`.

**Determinism.** No random number generation is used in estimation: A initializes from `spatialreg` fits, B constructs a deterministic multi-start set, and C starts from the constant vector 0.01 — so no seed is required and reruns are bit-reproducible up to BLAS/thread scheduling effects. Optimizer settings: B — `optim(method = "BFGS")` in tanh space, `reltol = 1e-8`, `maxit = 220` (retry 320); C — `optim(method = "L-BFGS-B")`, bounds ±0.99, `maxit = 500`; inner GLS: 50 iterations, tolerance 10⁻⁶ (both). Σ is symmetrized as (S + Sᵀ)/2 at every update; C additionally floors Σ's eigenvalues at 10⁻⁸ (`regularize_symmetric`), and B rejects any GLS iterate whose Σ has a non-positive eigenvalue.

**Parallel configuration on the benchmark machine.** A: PSOCK cluster, default `MSTR_CORES` = physical cores − 1 = 3 (overridable via env var `MSTR_CORES`), BLAS threads pinned to 1 per worker; B: per-model PSOCK, auto rule for n = 400 gives ⌊0.75 × (8 − 1)⌋ = 5 workers; C: `doParallel`, min(11, 8 − 1) = 7 workers.

**Input file integrity (SHA-256).**

```
7c0f63b34ab6d491241b7e277808e289f09f3d188026b0a94b523bc6eafc1fac  simulated_data_1111_n400_T5.csv
46269338737e6cbcd37ed60084b72793da6bbaee7538987b394f3390be4d235d  spatial_weights_n400.csv
98bb5b2c355849cecd5f9137adab25037700060da60195926aa1309d7dac464b  implement-A.r
764fdf056d89949b11099c7a9fa3248d99e9535881f3f48383e487fd80e48cfe  implement-B.r
dbbdfd44aa8e90d86170aee88b399558dcce74f792639a72b185d1c379a89ec7  implement-C.r
```

Checksums of the generated output files will be published with each tagged release.

**How to run.** All three scripts resolve the data files relative to the script directory as `../simulated_data_1111_n400_T5.csv` and `../spatial_weights_n400.csv`; place each script one directory below the data (or edit `DATA_FILE`/`data_file` at the top of the script).

```sh
# implement-B.r: non-interactive runs must set the output root policy
MSTR_USE_CURRENT_DIR=yes Rscript implement-B.r      # or MSTR_OUT_ROOT=/path/to/out

# implement-C.r
Rscript implement-C.r

# implement-A.r (requires spdep/spatialreg; worker count via MSTR_CORES)
MSTR_CORES=4 Rscript implement-A.r
```

**Criteria formulas as computed** (for this dataset): pAIC = −2ℓ(θ̂) + 2·deff, pBIC = −2ℓ(θ̂) + log(800)·deff, with deff of Section 2.4; ordinary AIC/BIC use the nominal parameter counts of the table in Section 5; the averaged pseudo-R̄² of Eq. (26) is reported descriptively by B.


## 7. Limitations and future work

(i) The runtime benchmark of Section 3.5 consists of a single end-to-end run of each script on one Windows workstation; replicated timings and per-stage profiling (optimization vs. Hessian evaluation vs. I/O) are needed before attributing the observed ratios to specific mechanisms. (ii) All statistical verification rests on a single simulated realization and on the conditional single-period estimation at t = 2; a Monte Carlo study (≥ 300 replications at n = 100/400/900) reporting true-model selection rates, over-/under-selection rates, parameter bias and RMSE, convergence rates, and the distribution of the selected γ is required before any claim about the selection performance of pAIC/pBIC, as is an extension to joint estimation over all T periods. (iii) The cross-implementation agreement of Section 3.3 covers four models within a common pipeline; a full end-to-end comparison of the three actual R scripts — all eleven models, the entire γ grid, identical inputs and initial values — remains to be run and published. (iv) The pBIC path of Section 3.4 terminated at the boundary of the examined grid; the grid must be extended (or the γ → ∞ limit compared directly) before the pBIC-selected γ is interpreted. (v) The Hessian-eigenvalue flooring and fallback rules of Section 2.4 are approximations; the per-fit diagnostics listed there (minimum eigenvalue, negative-eigenvalue count, condition number, deff before/after flooring, finite-difference step sensitivity) should be added to the standard output. (vi) Porting A's symmetric-polynomial log-determinant into B/C to lift the K = 2 restriction, and unifying the three programs' admissibility rules (Section 2.3), are natural next steps.
