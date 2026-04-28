# Methodology notes

This document supplements the main paper with implementation-level details
about the analytical pipeline. For the formal description, see the main paper
(`paper/main.tex`) Section 4 and the supplementary material
(`supplementary/supplementary.tex`) Section S6.

---

## 1. Avalanche identification

A stress event for country *c* at month *t* is defined by

```
Δlog s_{c,t} > σ × σ̂_c
```

where σ̂_c is the in-sample standard deviation of `Δlog s_{c,t}`. We use σ = 1.5
as baseline (T03 reports country thresholds in basis points). An avalanche of
size *S* at month *t* is the count of countries simultaneously stressed.

Threshold robustness over σ ∈ [1.25, 2.50] is documented in T08 and Figure 7.

## 2. Power-law estimation

Following Clauset, Shalizi and Newman (2009):

1. For each candidate `x_min ∈ {x_1, x_2, ..., x_n}`, compute the conditional
   MLE α̂(x_min).
2. Compute the Kolmogorov-Smirnov distance D between the empirical CCDF and
   the fitted power-law CCDF on the tail x ≥ x_min.
3. Select x̂_min minimizing D. Report the corresponding α̂.

Implemented via the `powerlaw` package (Alstott et al.).

## 3. Vuong test (PL vs LN, PL vs Exp)

For non-nested alternatives the normalized Vuong statistic is

```
R = Σ log[f1/f2] / (sqrt(n) × σ̂_R)
```

asymptotically standard normal under the null of equal fits. Following
Stumpf and Porter (2012), |R| < 1.5 indicates the heavy-tail boundary
regime.

## 4. Filtered networks

Three networks are computed on each rolling 24-month window:

| Filter | Definition | Edges |
|--------|-----------|-------|
| corr | |ρ_{ij}| > 0.4 | variable, up to n(n-1)/2 = 55 |
| partial | |ρ^{part}_{ij}| > 0.4 | variable |
| mst | minimum spanning tree on d_{ij} = 1 - |ρ_{ij}| | exactly n - 1 = 10 |

## 5. Network geometry

For each network and snapshot:

- Forman-Ricci curvature on each edge (Forman 2003 formula).
- Ollivier-Ricci curvature on each edge (Ollivier 2009): solves an exact
  earth-mover linear program of size |N(u)| × |N(v)| per edge with
  scipy.optimize.linprog (HiGHS solver).
- Spectral radius ρ in two normalizations: max-row-sum and Frobenius.
- Spectral Fragility Index: SFI = ρ / (1 - ρ).
- Network density and Herfindahl-Hirschman index of node strengths.

## 6. Placebo synchronization test

Reshuffling preserves country-level stress frequencies but breaks
cross-country synchronization. Specifically, for each bootstrap replication
*b*:

1. For each country *c*, independently permute the binary stress indicator
   `1[Δlog s_{c,t} > σ × σ̂_c]` across the time index.
2. Compute four synchronization statistics: max avalanche size, mean
   positive avalanche size, share P(S ≥ 4), share P(S ≥ 6).
3. The *p*-value is the proportion of placebo replications attaining or
   exceeding the observed statistic.

We use B = 1000 replications with seed 2026.

## 7. Lead-lag predictive analysis

For each network metric x_t at snapshot t:

1. Compute the maximum avalanche size in (t, t + 3]:

   ```
   y_t = max{S_u : u ∈ (t, t + 3]}
   ```

2. Pearson correlation r between x_t and y_t over all valid (t, y_t) pairs.

The 3-month horizon is set to match the rolling-window step of the network
construction.

## 8. Computational notes

- The Ollivier-Ricci LP is the bottleneck. Each snapshot has ~30-50 edges,
  each LP has up to 11 × 11 = 121 decision variables. Total runtime ~5-10
  minutes on a mid-range laptop.
- The placebo test is fully vectorizable but uses an explicit per-country
  permutation loop for clarity. ~30 seconds for B = 1000.
- All randomness uses np.random.default_rng(seed) for reproducibility.

## References

See `paper/references.bib` for the full bibliography. Key methodological
references:

- Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). *Power-law
  distributions in empirical data.* SIAM Review, 51(4), 661-703.
- Vuong, Q. H. (1989). *Likelihood ratio tests for model selection and
  non-nested hypotheses.* Econometrica, 57(2), 307-333.
- Stumpf, M. P. H., & Porter, M. A. (2012). *Critical truths about power
  laws.* Science, 335(6069), 665-666.
- Forman, R. (2003). *Bochner's method for cell complexes and combinatorial
  Ricci curvature.* Discrete & Computational Geometry, 29(3), 323-374.
- Ollivier, Y. (2009). *Ricci curvature of Markov chains on metric spaces.*
  Journal of Functional Analysis, 256(3), 810-864.
- Sandhu, R. S., Georgiou, T. T., & Tannenbaum, A. R. (2016). *Ricci
  curvature: An economic indicator for market fragility and systemic risk.*
  Science Advances, 2(5), e1501495.
- Mantegna, R. N. (1999). *Hierarchical structure in financial markets.*
  The European Physical Journal B, 11(1), 193-197.
