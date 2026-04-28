# Data dictionary

This document describes the columns in each CSV file produced by the
analytical pipeline (`code/01_pipeline_v3_real_embi.py`) and consumed by the
table generator (`code/02_generate_paper_tables.py`).

All files are in long format with one observation per row.

---

## EMBI_mensual.xlsx

The raw input. J.P. Morgan EMBI Global Diversified monthly stripped spreads,
sourced from a licensed Bloomberg terminal. Sheet `Hoja1`, header at row 7,
monthly block at columns 19-32. Spreads are in **percentage points**; the
pipeline multiplies by 100 to convert to **basis points**.

---

## latam_spread_panel_real_embi.csv

Cleaned long-format spread panel. *N* = 2,453 observations (11 countries ×
223 months).

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | Calendar month (first day of month) |
| country | str (3) | ISO-3 country code (ARG, BRA, CHL, COL, DOM, ECU, MEX, PAN, PER, SLV, URY) |
| spread_bp | float | EMBI Global Diversified stripped spread, basis points |

---

## latam_avalanche_events_real.csv

Country-stress events at the σ = 1.5 baseline threshold. *N* = 140 rows.

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | Calendar month of stress event |
| country | str (3) | ISO-3 country code |
| dlog_spread | float | Δlog s_{c,t} = log(s_{c,t}) - log(s_{c,t-1}) |
| overshoot | float | dlog_spread - σ × σ̂_c (overshoot above threshold) |
| jump_bp | float | Δs_{c,t} = s_{c,t} - s_{c,t-1}, basis points |
| spread_bp | float | s_{c,t} (level after the jump), basis points |

---

## soc_intereventimes_panel_real.csv

Country-level rolling-window estimates of the inter-event-time tail
exponent.

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | End of the rolling window |
| country | str (3) | ISO-3 country code |
| alpha_iet | float | Estimated tail exponent of inter-event times in the 60-month window |
| xmin_iet | float | Optimal x_min from KS minimization, in months |
| n_events_window | int | Number of stress events used in the window fit |
| mean_wait_months | float | Mean inter-event time, months |

Note: values may be NaN where fewer than 3 events fall in the window.

---

## fiscal_network_edges_real.csv

Edge list for all three filtered networks across all rolling-window
snapshots.

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | Snapshot date (end of 24-month window) |
| type | str | Network filter: `corr`, `partial`, or `mst` |
| src | str (3) | Source country (ISO-3) |
| dst | str (3) | Destination country (ISO-3) |
| weight | float | Edge weight \|ρ\| or \|ρ_partial\| |

---

## network_geometry_snapshots_{corr,partial,mst}_real.csv

Network-level summary metrics, one row per snapshot per filter type.

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | Snapshot date |
| type | str | `corr`, `partial`, or `mst` |
| kappa_ollivier | float | Mean Ollivier-Ricci curvature over edges |
| kappa_forman | float | Mean Forman-Ricci curvature over edges |
| rho_max | float | Spectral radius (max-row-sum normalization) |
| rho_frob | float | Spectral radius (Frobenius normalization) |
| sfi_max | float | Spectral Fragility Index using rho_max |
| sfi_frob | float | Spectral Fragility Index using rho_frob |
| hhi | float | Herfindahl-Hirschman index of node strengths |
| density | float | Network density |
| n_nodes | int | Number of nodes (always 11) |
| n_edges | int | Number of edges in the snapshot |

---

## network_geometry_country_{corr,partial,mst}_real.csv

Country-level summary metrics within each network snapshot.

| Column | Type | Description |
|--------|------|-------------|
| date | YYYY-MM-DD | Snapshot date |
| type | str | `corr`, `partial`, or `mst` |
| country | str (3) | ISO-3 country code |
| kappa_country | float | Mean Ollivier-Ricci curvature over edges incident on the country |
| strength | float | Sum of edge weights incident on the country |
| degree | int | Number of edges incident on the country in this snapshot |

---

## Country-code reference

| Code | Country | Notes |
|------|---------|-------|
| ARG | Argentina | Two sovereign credit events in sample (2014 holdout, 2020 default) |
| BRA | Brazil | Largest panel issuer; 2015-16 commodity stress, 2018 elections |
| CHL | Chile | Most stable; max spread ~390 bp |
| COL | Colombia | 2014-15 commodity stress, 2021 protests |
| DOM | Dominican Republic | Less liquid; lowest event count (n=7) |
| ECU | Ecuador | Defaults in 2008-09 and 2020 |
| MEX | Mexico | Investment-grade; deeply integrated with U.S. |
| PAN | Panama | Investment-grade; lost rating Mar 2024 |
| PER | Peru | Investment-grade; political volatility 2018-23 |
| SLV | El Salvador | 2017-18 pre-default stress, Bukele bond restructuring |
| URY | Uruguay | Most resilient post-2020; lowest spreads end-of-sample |
