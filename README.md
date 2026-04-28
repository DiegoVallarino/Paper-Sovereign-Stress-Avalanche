# latam-soc-embi

**Self-Organized Criticality in Latin American Sovereign Spreads:
Avalanche Dynamics, Network Geometry, and Fragility Regimes**

Replication package for an empirical analysis of self-organized criticality
(SOC) in Latin American sovereign credit risk using monthly J.P. Morgan EMBI
Global Diversified spreads (eleven LATAM economies, October 2007 to April
2026, *N* = 2,453 observations).

> Diego Vallarino, Ph.D.
> Counselor to the Board of Directors, Inter-American Development Bank and IDB Invest
> Sandpile Economics Program

---

## Abstract

**Purpose.** This paper tests for empirical signatures of self-organized
criticality (SOC) in Latin American sovereign credit risk and characterizes
the geometry of cross-country contagion using direct EMBI data.

**Design / methodology.** Monthly J.P. Morgan EMBI Global Diversified spreads
for eleven LATAM economies (October 2007 to April 2026, *N* = 2,453
observations) are used to build empirical avalanches via country-specific
stress thresholds. Power-law tail estimation follows Clauset, Shalizi and
Newman (2009) with Vuong tests against lognormal and exponential
alternatives. Three filtered networks (correlation, partial correlation,
minimum spanning tree) on rolling 24-month windows yield Forman-Ricci and
Ollivier-Ricci curvature plus a spectral fragility index (SFI). A *B* = 1,000
reshuffling placebo formally tests whether observed synchronization is
non-trivial.

**Findings.** Three independent statistical objects exhibit a heavy-tail
boundary regime in the sense of Stumpf and Porter (2012): avalanche size
(α̂ = 1.77), pooled |Δs| (α̂ = 2.37) and inter-event times (α̂ = 2.38).
Power-law cannot be rejected against lognormal at conventional significance
levels for any of the three objects, while it is preferred over an
exponential alternative. Country participation is heterogeneous and
dominated by mid-credibility, highly traded sovereigns (Brazil, Mexico,
Peru, Colombia, Panama), *not* by high-spread economies. The placebo
decisively rejects independence of country-level stress (*p* < 0.001).
Partial-correlation network spectral fragility predicts the maximum
3-month-ahead avalanche size (*r* = +0.28, *p* = 0.025).

**Originality / value.** This is the first empirical SOC analysis of LATAM
sovereign credit using EMBI data, with formal power-law tests, a
placebo-based identification of non-trivial synchronization, and a spectral
fragility indicator that is operational for surveillance. The reformulation
that criticality emerges from financially integrated mid-credibility
sovereigns---rather than from fiscally peripheral ones---reframes how regional
fragility is monitored.

**Keywords:** Self-organized criticality; sovereign spreads; EMBI; Latin
America; financial networks; Ricci curvature; spectral fragility; contagion.

**JEL codes:** G15, F34, C58, E44, F65.

---

## Repository structure

```
latam-soc-embi/
├── README.md                       This file
├── LICENSE                         MIT License
├── CITATION.cff                    Citation metadata
├── requirements.txt                Python dependencies
├── .gitignore
│
├── code/
│   ├── 01_pipeline_v3_real_embi.py     Main analytical pipeline
│   └── 02_generate_paper_tables.py     LaTeX tables generator
│
├── data/
│   ├── EMBI_mensual.xlsx                       Raw J.P. Morgan EMBI input **
│   ├── latam_spread_panel_real_embi.csv        Cleaned panel (long format)
│   ├── latam_avalanche_events_real.csv         Country-stress events
│   ├── soc_intereventimes_panel_real.csv       Inter-event-time analysis
│   ├── network_geometry_snapshots_corr_real.csv
│   ├── network_geometry_snapshots_partial_real.csv
│   ├── network_geometry_snapshots_mst_real.csv
│   ├── network_geometry_country_corr_real.csv
│   ├── network_geometry_country_partial_real.csv
│   ├── network_geometry_country_mst_real.csv
│   └── fiscal_network_edges_real.csv
│
├── figures/                         Figures in PNG format
│   ├── fig1_latam_spreads.png
│   ├── fig2_avalanche_soc.png
│   ├── fig3_network_snapshots.png
│   ├── fig4_geometry_dynamics.png
│   ├── fig5_country_stress_recurrence.png
│   ├── fig6_stress_heatmap.png
│   ├── fig7_threshold_robustness.png
│   └── fig8_dspread_diagnostic.png
│
├── tables/                          Twelve LaTeX tables (paper-ready)
│   ├── master.tex                       Master include
│   ├── T01_data_coverage.tex
│   ├── T02_avalanche_size_distribution.tex
│   ├── T03_country_stress_thresholds.tex
│   ├── T04_yearly_chronology.tex
│   ├── T05_top_avalanches.tex
│   ├── T06_top_individual_jumps.tex
│   ├── T07_tail_diagnostics.tex
│   ├── T08_threshold_robustness.tex
│   ├── T09_placebo_synchronization.tex
│   ├── T10_network_regime_comparison.tex
│   ├── T11_network_predictive_association.tex
│   └── T12_hypothesis_summary.tex
│
├── supplementary/                   Supplementary material LaTeX source
│   ├── supplementary.tex
│   ├── references.bib
│   ├── T*.tex
│   └── fig*.png
│
└── docs/
    ├── METHODOLOGY.md               Methodological notes
    └── DATA_DICTIONARY.md           CSV column documentation
** available upon reasonable request
```

## Quick start

### Software requirements
- Python 3.11 or later
- See `requirements.txt`. Install with:

```bash
pip install -r requirements.txt
```

The `powerlaw` package is required for the Clauset, Shalizi and Newman (2009)
maximum-likelihood power-law fit. If pip refuses on managed Python
installations, use:

```bash
pip install -r requirements.txt --break-system-packages
```

### Reproducing all results

The pipeline is split into two scripts. Each takes a few minutes on a
mid-range laptop.

```bash
# 1. Run the analytical pipeline (loads EMBI data, builds networks,
#    computes geometry, generates 8 figures and 10 CSV outputs).
#    Edit BASE_DIR and EMBI_PATH at the top of the script first.
python code/01_pipeline_v3_real_embi.py

# 2. Generate the 12 paper-ready LaTeX tables.
#    Edit DATA_DIR and OUT_DIR at the top of the script first.
python code/02_generate_paper_tables.py
```

### Compiling the paper

```bash
cd paper/
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

### Compiling the supplementary material

```bash
cd supplementary/
pdflatex supplementary.tex
bibtex supplementary
pdflatex supplementary.tex
pdflatex supplementary.tex
```

## Hypotheses tested in the paper

| #   | Hypothesis                                                                                       | Outcome              |
|-----|--------------------------------------------------------------------------------------------------|----------------------|
| H1  | Country-level participation in stress avalanches is heterogeneous.                              | Confirmed            |
| H2  | Mid-credibility, highly traded sovereigns dominate avalanche participation.                     | Confirmed            |
| H3  | Avalanche size, |Δs|, and inter-event times all exhibit heavy-tail boundary regime.            | Confirmed            |
| H4  | Cross-country stress synchronization is non-trivial.                                             | Confirmed (*p*<0.001)|
| H5  | Network geometry shifts around large stress months and predicts near-future avalanche size.     | Partially confirmed  |
| H6  | Threshold robustness: avalanche evidence is robust to choice of σ.                              | Confirmed            |

## Key empirical findings

- The avalanche-size distribution exhibits a power-law tail with α̂ = 1.77,
  inside the canonical SOC range [1, 3].
- Three universal-stress events (S = 11) coincide exactly with the Lehman
  collapse (September 2008), the post-Lehman EM rout (October 2008), and the
  COVID-19 global rout (March 2020).
- The country-event ranking is dominated by Brazil (18 events), Mexico (16),
  Peru (16), Colombia (14) and Panama (14). Argentina and Ecuador have only
  10 events each despite the highest spread levels.
- The reshuffling placebo decisively rejects independence: the observed
  *S*<sub>max</sub> = 11 is unattainable in 1,000 placebo replications.
- Partial-correlation Frobenius spectral fragility correlates with maximum
  3-month-ahead avalanche size at *r* = +0.28 (*p* = 0.025).

## Data sources

The primary data source is the **J.P. Morgan EMBI Global Diversified** monthly
stripped spread series, obtained from a licensed Bloomberg terminal. The raw
file is included as `data/EMBI_mensual.xlsx` for reproducibility.

## Citation

If you use this code or data, please cite:

```bibtex
@unpublished{Vallarino2026SOC,
  author = {Vallarino, Diego},
  title  = {Self-Organized Criticality in {L}atin {A}merican Sovereign Spreads:
            Avalanche Dynamics, Network Geometry, and Fragility Regimes},
  year   = {2026},
  note   = {Working paper},
  url    = {https://github.com/diegovallarino/latam-soc-embi}
}
```

See also `CITATION.cff` for GitHub-friendly citation metadata.

## License

This repository is released under the MIT License (see `LICENSE`).
The underlying J.P. Morgan EMBI data are reproduced for academic
replication purposes only and remain the property of J.P. Morgan.

## Contact

Diego Vallarino, Ph.D.
Counselor to the Board of Directors
Inter-American Development Bank and IDB Invest
Email: diegoval@iadb.org

The views expressed in this work are those of the author and do not
necessarily reflect those of the Inter-American Development Bank or its
Board of Directors.
