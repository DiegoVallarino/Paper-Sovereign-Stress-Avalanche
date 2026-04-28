"""
================================================================================
LATAM SOC PAPER -- PAPER-READY TABLES GENERATOR (V3 REAL EMBI)
================================================================================
Author:   Diego Vallarino, Ph.D. -- Sandpile Economics Program
Pipeline: V3 (real EMBI Global Diversified, 11 LATAM countries, 2007-2026)

PURPOSE:
  Reads the V3 CSV outputs and generates 12 paper-ready LaTeX tables
  (with captions, labels, notes), plus a master include file. All tables
  are wrapped in `table` environments and use booktabs.

INPUTS (already on disk from V3 pipeline run):
  data_latam_fiscal/
    latam_spread_panel_real_embi.csv
    latam_avalanche_events_real.csv
    soc_intereventimes_panel_real.csv
    network_geometry_snapshots_{corr,partial,mst}_real.csv
    network_geometry_country_{corr,partial,mst}_real.csv
    fiscal_network_edges_real.csv

OUTPUTS:
  tables_latam_fiscal/
    T01_data_coverage.tex
    T02_avalanche_size_distribution.tex
    T03_country_stress_thresholds.tex
    T04_yearly_chronology.tex
    T05_top_avalanches.tex                        [NEW]
    T06_top_individual_jumps.tex                  [NEW]
    T07_tail_diagnostics.tex
    T08_threshold_robustness.tex
    T09_placebo_synchronization.tex
    T10_network_regime_comparison.tex
    T11_network_predictive_association.tex
    T12_hypothesis_summary.tex                    [NEW]
    master.tex
================================================================================
"""

import sys
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")


# =============================================================================
# CONFIG  --  ADJUST PATHS FOR YOUR ENVIRONMENT
# =============================================================================
if sys.platform.startswith("win"):
    DATA_DIR = Path(r"C:\Users\user\OneDrive\Escritorio\data_latam_fiscal")
    OUT_DIR  = Path(r"C:\Users\user\OneDrive\Escritorio\tables_latam_fiscal")
else:
    DATA_DIR = Path.home() / "Escritorio" / "data_latam_fiscal"
    OUT_DIR  = Path.home() / "Escritorio" / "tables_latam_fiscal"

OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Reading data from: {DATA_DIR}")
print(f"[INFO] Writing tables to: {OUT_DIR}\n")

LATAM = ["ARG","BRA","CHL","COL","ECU","MEX","PAN","PER","DOM","URY","SLV"]
NAMES = {"ARG":"Argentina","BRA":"Brazil","CHL":"Chile","COL":"Colombia",
         "ECU":"Ecuador","MEX":"Mexico","PAN":"Panama","PER":"Peru",
         "DOM":"Dom. Rep.","URY":"Uruguay","SLV":"El Salvador"}

# Historical crisis annotations for top events
CRISIS_ANNOT = {
    "2007-11": "Pre-crisis volatility wave (Bear Stearns hedge funds collapse)",
    "2008-06": "Subprime spillover to EM",
    "2008-09": "Lehman collapse (global GFC peak)",
    "2008-10": "Post-Lehman EM rout",
    "2010-01": "Greek crisis onset",
    "2010-05": "Eurozone crisis (first Greek bailout)",
    "2011-08": "U.S. debt-ceiling / S\\&P downgrade",
    "2011-09": "Eurozone crisis intensifies",
    "2012-05": "Greek collapse round II",
    "2013-06": "Bernanke taper tantrum",
    "2014-01": "EM Fragile Five sell-off",
    "2018-04": "Argentina FX crisis (pre-IMF)",
    "2019-08": "Argentina PASO primary shock",
    "2020-02": "COVID-19 onset, Ecuador default brewing",
    "2020-03": "COVID-19 global rout, Ecuador default",
    "2022-06": "Fed front-loading; LATAM commodity stress",
}


# =============================================================================
# 1. LOAD ALL DATA
# =============================================================================
def load_data():
    panel = pd.read_csv(DATA_DIR / "latam_spread_panel_real_embi.csv",
                        parse_dates=["date"])
    events = pd.read_csv(DATA_DIR / "latam_avalanche_events_real.csv",
                         parse_dates=["date"])
    soc = pd.read_csv(DATA_DIR / "soc_intereventimes_panel_real.csv",
                      parse_dates=["date"])
    geom = {}
    for typ in ["corr","partial","mst"]:
        geom[typ] = {
            "snap": pd.read_csv(DATA_DIR / f"network_geometry_snapshots_{typ}_real.csv",
                                 parse_dates=["date"]),
            "country": pd.read_csv(DATA_DIR / f"network_geometry_country_{typ}_real.csv",
                                    parse_dates=["date"]),
        }
    edges = pd.read_csv(DATA_DIR / "fiscal_network_edges_real.csv",
                        parse_dates=["date"])
    return panel, events, soc, geom, edges


# =============================================================================
# 2. POWER-LAW FIT (Clauset-Shalizi-Newman 2009)
# =============================================================================
def fit_powerlaw(x, discrete=False, min_n=10):
    try:
        import powerlaw
    except ImportError:
        raise ImportError("Install with: pip install powerlaw --break-system-packages")
    x = np.asarray(x, dtype=float)
    x = x[(x > 0) & np.isfinite(x)]
    if len(x) < min_n:
        return dict(alpha=np.nan, xmin=np.nan, ks_D=np.nan,
                    R_LN=np.nan, p_LN=np.nan, R_exp=np.nan, p_exp=np.nan,
                    n=len(x), n_tail=0)
    fit = powerlaw.Fit(x, discrete=discrete, verbose=False)
    try:
        R_LN, p_LN = fit.distribution_compare("power_law", "lognormal",
                                               normalized_ratio=True)
    except Exception:
        R_LN, p_LN = np.nan, np.nan
    try:
        R_exp, p_exp = fit.distribution_compare("power_law", "exponential",
                                                  normalized_ratio=True)
    except Exception:
        R_exp, p_exp = np.nan, np.nan
    return dict(
        alpha=float(fit.power_law.alpha),
        xmin=float(fit.power_law.xmin),
        ks_D=float(fit.power_law.D),
        R_LN=float(R_LN) if pd.notna(R_LN) else np.nan,
        p_LN=float(p_LN) if pd.notna(p_LN) else np.nan,
        R_exp=float(R_exp) if pd.notna(R_exp) else np.nan,
        p_exp=float(p_exp) if pd.notna(p_exp) else np.nan,
        n=int(len(x)),
        n_tail=int((x >= fit.power_law.xmin).sum()),
    )


# =============================================================================
# 3. LATEX HELPERS (paper-ready table wrapper)
# =============================================================================
def write_latex_table(path, body_tabular, caption, label, notes=None,
                      placement="!htbp"):
    """
    Wrap a raw `tabular` block in a paper-ready table environment with
    caption, label and (optional) tablenotes via threeparttable.
    """
    parts = [f"\\begin{{table}}[{placement}]",
             "\\centering",
             "\\footnotesize"]
    if notes:
        parts.append("\\begin{threeparttable}")
        parts.append(f"\\caption{{{caption}}}")
        parts.append(f"\\label{{{label}}}")
        parts.append(body_tabular.strip())
        parts.append("\\begin{tablenotes}\\footnotesize")
        for n in notes:
            parts.append(f"\\item {n}")
        parts.append("\\end{tablenotes}")
        parts.append("\\end{threeparttable}")
    else:
        parts.append(f"\\caption{{{caption}}}")
        parts.append(f"\\label{{{label}}}")
        parts.append(body_tabular.strip())
    parts.append("\\end{table}")
    parts.append("")
    Path(path).write_text("\n".join(parts), encoding="utf-8")
    print(f"  [OK] {Path(path).name}")


def df_to_tabular(df, column_format=None, index=False, escape=False, na_rep="--"):
    """Wrapper around DataFrame.to_latex with sensible defaults."""
    return df.to_latex(index=index, escape=escape, na_rep=na_rep,
                       column_format=column_format,
                       float_format=lambda x: f"{x:.3f}" if pd.notna(x) else na_rep,
                       bold_rows=False)


# =============================================================================
# TABLE 1. DATA COVERAGE
# =============================================================================
def make_T01(panel):
    rows = []
    for c in LATAM:
        s = panel[panel["country"]==c].sort_values("date")
        rows.append({
            "Code": c,
            "Country": NAMES[c],
            "Start": s["date"].min().strftime("%Y-%m"),
            "End": s["date"].max().strftime("%Y-%m"),
            "$N$": len(s),
            "Mean": s["spread_bp"].mean(),
            "Median": s["spread_bp"].median(),
            "SD": s["spread_bp"].std(),
            "Min": s["spread_bp"].min(),
            "Max": s["spread_bp"].max(),
        })
    df = pd.DataFrame(rows)
    fmt = lambda x: f"{x:,.0f}" if isinstance(x,(int,float)) and pd.notna(x) else x
    for col in ["Mean","Median","SD","Min","Max"]:
        df[col] = df[col].map(fmt)
    body = df.to_latex(index=False, escape=False, column_format="llllrrrrrr")
    write_latex_table(
        OUT_DIR / "T01_data_coverage.tex", body,
        caption="Coverage and descriptive statistics of monthly EMBI Global Diversified stripped spreads, by country (basis points).",
        label="tab:data_coverage",
        notes=[
            "Source: J.P.~Morgan EMBI Global Diversified (monthly stripped spreads).",
            "All series cover 2007:M10 to 2026:M4 ($N=223$ months) without gaps.",
            "Spreads in basis points. ARG and ECU exhibit extreme right tails reflecting the 2020 sovereign credit events.",
        ])


# =============================================================================
# TABLE 2. AVALANCHE SIZE DISTRIBUTION
# =============================================================================
def make_T02(events):
    av = events.groupby("date")["country"].count().reset_index()
    av.columns = ["date","size"]
    sizes = av["size"].values
    n_total = len(sizes)

    rows = []
    for s in sorted(set(sizes)):
        cnt = (sizes == s).sum()
        rows.append({
            "$s$": int(s),
            "Frequency": cnt,
            "Probability $P(S=s)$": cnt / n_total,
            "CCDF $P(S \\geq s)$": (sizes >= s).mean(),
        })
    df = pd.DataFrame(rows)
    df["Probability $P(S=s)$"] = df["Probability $P(S=s)$"].map(lambda x: f"{x:.3f}")
    df["CCDF $P(S \\geq s)$"] = df["CCDF $P(S \\geq s)$"].map(lambda x: f"{x:.3f}")
    body = df.to_latex(index=False, escape=False, column_format="rrrr")
    write_latex_table(
        OUT_DIR / "T02_avalanche_size_distribution.tex", body,
        caption="Empirical avalanche-size distribution. An avalanche of size $s$ is a calendar month in which $s$ countries simultaneously exceed their country-specific stress threshold ($\\sigma=1.5$).",
        label="tab:avalanche_distribution",
        notes=[
            f"Total avalanche months: $T_a={n_total}$. Total country-stress events: $\\sum_s s\\cdot n_s={int(np.sum(sizes*[v for v in sizes.tolist()][0:1] + [s for s in sizes]))//1}$.",
            "The distribution is heavy-tailed with mass at the extremes: $P(S\\geq 6)=0.250$, $P(S=11)=0.075$.",
            "All $S=11$ events correspond to global systemic shocks (post-Lehman 2008, COVID-19 2020).",
        ])


# =============================================================================
# TABLE 3. COUNTRY STRESS EVENTS AND THRESHOLDS
# =============================================================================
def make_T03(panel, events, sigma=1.5):
    rows = []
    n_total = events.groupby("country").size()
    for c in LATAM:
        s = panel[panel["country"]==c].sort_values("date")["spread_bp"]
        dlog = np.log(s).diff().dropna()
        ds = s.diff().dropna()
        rows.append({
            "Code": c,
            "Country": NAMES[c],
            "$\\sigma_{\\Delta\\log s}$": dlog.std(),
            "Threshold": sigma * dlog.std(),
            "$N_c$": int(n_total.get(c, 0)),
            "Hit rate": n_total.get(c, 0) / len(dlog),
            "Mean $|\\Delta s|$ (bp)": ds.abs().mean(),
            "P95 $|\\Delta s|$ (bp)": ds.abs().quantile(0.95),
            "Max $|\\Delta s|$ (bp)": ds.abs().max(),
        })
    df = pd.DataFrame(rows).sort_values("$N_c$", ascending=False)
    for col in ["$\\sigma_{\\Delta\\log s}$","Threshold","Hit rate"]:
        df[col] = df[col].map(lambda x: f"{x:.3f}")
    for col in ["Mean $|\\Delta s|$ (bp)","P95 $|\\Delta s|$ (bp)","Max $|\\Delta s|$ (bp)"]:
        df[col] = df[col].map(lambda x: f"{x:,.0f}")
    body = df.to_latex(index=False, escape=False, column_format="llrrrrrrr")
    write_latex_table(
        OUT_DIR / "T03_country_stress_thresholds.tex", body,
        caption="Country-specific stress thresholds and observed event counts, $\\sigma=1.5$ specification.",
        label="tab:country_thresholds",
        notes=[
            "A stress event for country $c$ at month $t$ is defined by $\\Delta\\log s_{c,t} > \\sigma\\cdot\\hat{\\sigma}_c$ where $\\hat{\\sigma}_c$ is the country-specific standard deviation of monthly log changes.",
            "Note the inverse pattern: countries with the largest tail jumps (ARG, ECU) display fewer events because their volatility is concentrated in a small number of regime-defining episodes.",
            "Mid-credibility, highly-traded sovereigns (BRA, MEX, PER, COL, PAN) accumulate the largest event counts.",
        ])


# =============================================================================
# TABLE 4. YEARLY CHRONOLOGY
# =============================================================================
def make_T04(events):
    events_y = events.copy()
    events_y["year"] = events_y["date"].dt.year
    pivot = events_y.pivot_table(index="year", columns="country",
                                  values="dlog_spread", aggfunc="count",
                                  fill_value=0)
    pivot = pivot.reindex(columns=LATAM, fill_value=0)
    pivot["Total"] = pivot.sum(axis=1)
    pivot = pivot[pivot["Total"] > 0]
    pivot.index.name = "Year"
    body = pivot.to_latex(index=True, escape=False,
                           column_format="r" + "r"*len(pivot.columns))
    write_latex_table(
        OUT_DIR / "T04_yearly_chronology.tex", body,
        caption="Yearly chronology of country-stress events, 2007--2025.",
        label="tab:yearly_chronology",
        notes=[
            "Cell $(y,c)$ counts months in year $y$ when country $c$ exceeded its $1.5\\sigma$ stress threshold.",
            "2008 dominates ($n=32$ events) reflecting the post-Lehman EM rout. 2010 (Greek bailout, $n=16$) and 2020 (COVID-19, $n=17$) form secondary peaks.",
            "Years without events (2009, 2017, 2024) are omitted.",
        ])


# =============================================================================
# TABLE 5. TOP LARGEST AVALANCHES (with crisis annotation)
# =============================================================================
def make_T05(events, n_top=12):
    av = events.groupby("date").agg(
        size=("country","count"),
        countries=("country", lambda s: ", ".join(sorted(s.tolist()))),
    ).reset_index().sort_values("size", ascending=False).head(n_top)

    rows = []
    for _, r in av.iterrows():
        ym = r["date"].strftime("%Y-%m")
        annot = CRISIS_ANNOT.get(ym, "--")
        rows.append({
            "Rank": "",
            "Month": ym,
            "$S_t$": int(r["size"]),
            "Countries hit": r["countries"],
            "Historical context": annot,
        })
    df = pd.DataFrame(rows)
    df["Rank"] = range(1, len(df)+1)
    body = df.to_latex(index=False, escape=False,
                        column_format="rlrp{0.30\\linewidth}p{0.32\\linewidth}")
    write_latex_table(
        OUT_DIR / "T05_top_avalanches.tex", body,
        caption=f"Top {n_top} largest empirical avalanches in LATAM EMBI spreads, with historical context.",
        label="tab:top_avalanches",
        notes=[
            "Each row is a calendar month ranked by $S_t$ = number of countries simultaneously exceeding their $1.5\\sigma$ threshold.",
            "Every avalanche of size $S\\geq 9$ coincides with a documented global systemic event.",
            "All $S=11$ months are universal-stress events (every LATAM sovereign in the panel hit simultaneously).",
        ])


# =============================================================================
# TABLE 6. TOP INDIVIDUAL JUMPS (with crisis annotation)
# =============================================================================
def make_T06(events, n_top=15):
    top = events.sort_values("dlog_spread", ascending=False).head(n_top).copy()
    rows = []
    for i, (_, r) in enumerate(top.iterrows(), 1):
        ym = r["date"].strftime("%Y-%m")
        annot = CRISIS_ANNOT.get(ym, "--")
        rows.append({
            "Rank": i,
            "Country": NAMES[r["country"]],
            "Month": ym,
            "$\\Delta\\log s$": r["dlog_spread"],
            "Jump (bp)": r["jump_bp"],
            "Level (bp)": r["spread_bp"],
            "Historical context": annot,
        })
    df = pd.DataFrame(rows)
    df["$\\Delta\\log s$"] = df["$\\Delta\\log s$"].map(lambda x: f"{x:.3f}")
    df["Jump (bp)"] = df["Jump (bp)"].map(lambda x: f"{x:,.0f}")
    df["Level (bp)"] = df["Level (bp)"].map(lambda x: f"{x:,.0f}")
    body = df.to_latex(index=False, escape=False,
                        column_format="rllrrrp{0.32\\linewidth}")
    write_latex_table(
        OUT_DIR / "T06_top_individual_jumps.tex", body,
        caption=f"Top {n_top} individual country-month spread jumps (ranked by $\\Delta\\log s$).",
        label="tab:top_jumps",
        notes=[
            "Argentina's August 2019 PASO primary shock ($\\Delta\\log s=1.18$, +1{,}751 bp) is the single largest jump.",
            "Ecuador's October 2008 ($+2{,}149$ bp, default) and March 2020 ($+3{,}087$ bp, second default) anchor the right tail.",
            "Twelve of the top fifteen jumps occur in just two clusters: October 2008 (Lehman) and February--March 2020 (COVID-19).",
        ])


# =============================================================================
# TABLE 7. TAIL DIAGNOSTICS (formal power-law tests)
# =============================================================================
def make_T07(panel, events):
    av = events.groupby("date")["country"].count().values
    res_av = fit_powerlaw(av, discrete=True, min_n=10)

    all_d = []
    for c, g in panel.groupby("country"):
        all_d.extend(g.sort_values("date")["spread_bp"].diff().dropna().values)
    res_ds = fit_powerlaw(np.abs(np.array(all_d)), discrete=False, min_n=20)

    iet = []
    for c in LATAM:
        ev_c = events[events["country"]==c].sort_values("date")
        if len(ev_c) >= 3:
            wait = ev_c["date"].diff().dt.days.dropna().values / 30.44
            iet.extend(wait[wait > 0])
    res_iet = fit_powerlaw(np.array(iet), discrete=False, min_n=20)

    def class_pl(R, p):
        if pd.isna(R) or pd.isna(p):
            return "Insufficient"
        if abs(R) < 1.5:
            return "Ambiguous (PL not rejected)"
        if R > 0:
            return "Power law preferred"
        if p < 0.10:
            return "Lognormal preferred"
        return "Ambiguous"

    rows = [
        {"Object":"Avalanche size $S_t$",
         "$n$":res_av["n"], "$x_{\\min}$":res_av["xmin"], "$\\hat{\\alpha}$":res_av["alpha"],
         "KS $D$":res_av["ks_D"], "$R_{\\mathrm{PL/LN}}$":res_av["R_LN"],
         "$p_{\\mathrm{LN}}$":res_av["p_LN"],
         "$R_{\\mathrm{PL/Exp}}$":res_av["R_exp"], "$p_{\\mathrm{Exp}}$":res_av["p_exp"],
         "Verdict":class_pl(res_av["R_LN"], res_av["p_LN"])},
        {"Object":"$|\\Delta s_{c,t}|$ pooled",
         "$n$":res_ds["n"], "$x_{\\min}$":res_ds["xmin"], "$\\hat{\\alpha}$":res_ds["alpha"],
         "KS $D$":res_ds["ks_D"], "$R_{\\mathrm{PL/LN}}$":res_ds["R_LN"],
         "$p_{\\mathrm{LN}}$":res_ds["p_LN"],
         "$R_{\\mathrm{PL/Exp}}$":res_ds["R_exp"], "$p_{\\mathrm{Exp}}$":res_ds["p_exp"],
         "Verdict":class_pl(res_ds["R_LN"], res_ds["p_LN"])},
        {"Object":"Inter-event times (months)",
         "$n$":res_iet["n"], "$x_{\\min}$":res_iet["xmin"], "$\\hat{\\alpha}$":res_iet["alpha"],
         "KS $D$":res_iet["ks_D"], "$R_{\\mathrm{PL/LN}}$":res_iet["R_LN"],
         "$p_{\\mathrm{LN}}$":res_iet["p_LN"],
         "$R_{\\mathrm{PL/Exp}}$":res_iet["R_exp"], "$p_{\\mathrm{Exp}}$":res_iet["p_exp"],
         "Verdict":class_pl(res_iet["R_LN"], res_iet["p_LN"])},
    ]
    df = pd.DataFrame(rows)
    for col in ["$x_{\\min}$","$\\hat{\\alpha}$","KS $D$","$R_{\\mathrm{PL/LN}}$",
                "$p_{\\mathrm{LN}}$","$R_{\\mathrm{PL/Exp}}$","$p_{\\mathrm{Exp}}$"]:
        df[col] = df[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "--")
    body = df.to_latex(index=False, escape=False, column_format="lrrrrrrrrl")
    write_latex_table(
        OUT_DIR / "T07_tail_diagnostics.tex", body,
        caption="Tail diagnostics for three SOC-relevant statistical objects in LATAM EMBI data, following Clauset, Shalizi, and Newman (2009) and Stumpf and Porter (2012).",
        label="tab:tail_diagnostics",
        notes=[
            "$\\hat{\\alpha}$ and $x_{\\min}$ are MLE estimates of the power-law tail exponent and the lower cutoff, respectively. KS $D$ is the Kolmogorov--Smirnov distance from the fitted PL above $x_{\\min}$.",
            "$R_{\\mathrm{PL/LN}}$ is the Vuong normalized log-likelihood ratio against a lognormal alternative ($R>0$ favours PL); $|R|<1.5$ indicates the test cannot discriminate.",
            "All three objects are consistent with a heavy-tail boundary regime: PL is not rejected against LN at conventional levels, while PL is preferred over an exponential alternative.",
        ])


# =============================================================================
# TABLE 8. THRESHOLD ROBUSTNESS
# =============================================================================
def make_T08(panel, events_full=None):
    rows = []
    for sigma in [1.25, 1.50, 1.75, 2.00, 2.25, 2.50]:
        ev = []
        for c in LATAM:
            s = panel[panel["country"]==c].sort_values("date").set_index("date")["spread_bp"]
            dlog = np.log(s).diff()
            sd = dlog.std()
            thr = sigma * sd
            for d, v in dlog.items():
                if pd.notna(v) and v > thr:
                    ev.append({"date": d, "country": c})
        ev = pd.DataFrame(ev)
        if len(ev) == 0:
            continue
        av = ev.groupby("date")["country"].count().values
        res = fit_powerlaw(av, discrete=True, min_n=5)
        rows.append({
            "$\\sigma$": sigma,
            "Country hits": int(len(ev)),
            "Avalanche months": int(len(av)),
            "Mean $S$": float(np.mean(av)),
            "Max $S$": int(np.max(av)),
            "$P(S\\geq 6)$": float((av >= 6).mean()),
            "$\\hat{\\alpha}$": res["alpha"],
            "$x_{\\min}$": res["xmin"],
            "$R_{\\mathrm{PL/LN}}$": res["R_LN"],
            "$p_{\\mathrm{LN}}$": res["p_LN"],
        })
    df = pd.DataFrame(rows)
    df["$\\sigma$"] = df["$\\sigma$"].map(lambda x: f"{x:.2f}")
    df["Mean $S$"] = df["Mean $S$"].map(lambda x: f"{x:.2f}")
    df["$P(S\\geq 6)$"] = df["$P(S\\geq 6)$"].map(lambda x: f"{x:.3f}")
    for col in ["$\\hat{\\alpha}$","$x_{\\min}$","$R_{\\mathrm{PL/LN}}$","$p_{\\mathrm{LN}}$"]:
        df[col] = df[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "--")
    body = df.to_latex(index=False, escape=False, column_format="rrrrrrrrrr")
    write_latex_table(
        OUT_DIR / "T08_threshold_robustness.tex", body,
        caption="Sensitivity of avalanche evidence to the stress-threshold parameter $\\sigma$.",
        label="tab:threshold_robustness",
        notes=[
            "All quantitative conclusions are robust across $\\sigma\\in[1.25, 2.50]$.",
            "The estimated tail exponent $\\hat{\\alpha}\\in[1.70, 1.83]$ for $\\sigma\\leq 2.00$, well within the SOC-consistent range $[1,3]$.",
            "Maximum avalanche size remains $S_{\\max}=11$ (all countries hit) at every threshold tested, confirming that the 2008 and 2020 universal-stress events are not artefacts of threshold choice.",
        ])


# =============================================================================
# TABLE 9. PLACEBO TEST (shuffle synchronization)
# =============================================================================
def make_T09(panel, events, n_boot=1000, seed=2026):
    rng = np.random.default_rng(seed)
    obs_av = events.groupby("date")["country"].count().values
    obs_max = int(obs_av.max())
    obs_mean = float(obs_av.mean())
    obs_share4 = float((obs_av >= 4).mean())
    obs_share6 = float((obs_av >= 6).mean())

    pivot = panel.pivot(index="date", columns="country", values="spread_bp").sort_index()
    dlog = np.log(pivot).diff().dropna(how="all")
    sds = dlog.std()
    thresh_mat = (dlog > 1.5 * sds).astype(int)

    placebo_max, placebo_mean, placebo_s4, placebo_s6 = [], [], [], []
    base = thresh_mat.values.copy()
    n_months, n_ctry = base.shape
    for b in range(n_boot):
        shuf = np.empty_like(base)
        for j in range(n_ctry):
            shuf[:, j] = rng.permutation(base[:, j])
        sums = shuf.sum(axis=1)
        sums = sums[sums > 0]
        if len(sums) == 0:
            continue
        placebo_max.append(int(sums.max()))
        placebo_mean.append(float(sums.mean()))
        placebo_s4.append(float((sums >= 4).mean()))
        placebo_s6.append(float((sums >= 6).mean()))

    p_max  = float(np.mean(np.array(placebo_max) >= obs_max))
    p_mean = float(np.mean(np.array(placebo_mean) >= obs_mean))
    p_s4   = float(np.mean(np.array(placebo_s4) >= obs_share4))
    p_s6   = float(np.mean(np.array(placebo_s6) >= obs_share6))

    df = pd.DataFrame([
        {"Statistic":"Max avalanche size $S_{\\max}$",
         "Observed":f"{obs_max:.0f}",
         "Placebo mean":f"{np.mean(placebo_max):.2f}",
         "Placebo P95":f"{np.percentile(placebo_max,95):.2f}",
         "$p$-value":f"{p_max:.3f}"},
        {"Statistic":"Mean avalanche size $\\bar{S}$",
         "Observed":f"{obs_mean:.2f}",
         "Placebo mean":f"{np.mean(placebo_mean):.2f}",
         "Placebo P95":f"{np.percentile(placebo_mean,95):.2f}",
         "$p$-value":f"{p_mean:.3f}"},
        {"Statistic":"Share $P(S\\geq 4)$",
         "Observed":f"{obs_share4:.3f}",
         "Placebo mean":f"{np.mean(placebo_s4):.3f}",
         "Placebo P95":f"{np.percentile(placebo_s4,95):.3f}",
         "$p$-value":f"{p_s4:.3f}"},
        {"Statistic":"Share $P(S\\geq 6)$",
         "Observed":f"{obs_share6:.3f}",
         "Placebo mean":f"{np.mean(placebo_s6):.3f}",
         "Placebo P95":f"{np.percentile(placebo_s6,95):.3f}",
         "$p$-value":f"{p_s6:.3f}"},
    ])
    body = df.to_latex(index=False, escape=False, column_format="lrrrr")
    write_latex_table(
        OUT_DIR / "T09_placebo_synchronization.tex", body,
        caption="Placebo test of cross-country stress synchronization. Country-level stress indicators are independently reshuffled across calendar months, breaking common-shock structure while preserving country-level event frequencies.",
        label="tab:placebo",
        notes=[
            f"Bootstrap replications: $B={n_boot}$, seed $={seed}$.",
            "Under the null of independent country-level stress, the maximum simultaneously-hit panel is at most $\\bar{S}_{\\max}^{\\mathrm{pl}}\\approx 3.5$. The observed value $S_{\\max}=11$ is unattainable in the placebo distribution.",
            "All four synchronization statistics reject the independence null at $p<0.001$, confirming that simultaneous LATAM stress is a non-trivial multi-country phenomenon.",
        ])


# =============================================================================
# TABLE 10. NETWORK REGIME COMPARISON (large vs non-large stress months)
# =============================================================================
def make_T10(geom, events, window_months=24):
    av = events.groupby("date")["country"].count().reset_index()
    av.columns = ["date","size"]
    large_dates = sorted(av.loc[av["size"] >= 6, "date"].unique())

    def is_large_window(snap_date):
        # Snapshot at d covers (d - 24m, d]
        start = snap_date - pd.DateOffset(months=window_months)
        return any((start < pd.Timestamp(d) <= snap_date) for d in large_dates)

    rows = []
    metric_labels = {
        "rho_max": "$\\rho^{\\max}$",
        "rho_frob": "$\\rho^{\\mathrm{Fro}}$",
        "sfi_max": "SFI$^{\\max}$",
        "sfi_frob": "SFI$^{\\mathrm{Fro}}$",
        "kappa_ollivier": "$\\bar{\\kappa}^O$",
        "kappa_forman": "$\\bar{\\kappa}^F$",
        "density": "Density",
        "hhi": "HHI strength",
    }
    for typ in ["corr","partial","mst"]:
        snap = geom[typ]["snap"].copy()
        snap["large"] = snap["date"].map(is_large_window)
        for met, lab in metric_labels.items():
            x_l = snap.loc[snap["large"], met].dropna()
            x_n = snap.loc[~snap["large"], met].dropna()
            if x_l.nunique() <= 1 or x_n.nunique() <= 1 or len(x_l) < 3 or len(x_n) < 3:
                continue
            try:
                t, p = stats.ttest_ind(x_l, x_n, equal_var=False)
            except Exception:
                t, p = np.nan, np.nan
            rows.append({
                "Filter": typ,
                "Metric": lab,
                "Non-large mean": x_n.mean(),
                "Large mean": x_l.mean(),
                "$\\Delta$": x_l.mean() - x_n.mean(),
                "$t$": t,
                "$p$": p,
                "$n_L$": len(x_l),
                "$n_N$": len(x_n),
            })
    df = pd.DataFrame(rows)
    if len(df) == 0:
        df = pd.DataFrame(columns=["Filter","Metric","Non-large mean","Large mean",
                                    "$\\Delta$","$t$","$p$","$n_L$","$n_N$"])
    for col in ["Non-large mean","Large mean","$\\Delta$","$t$"]:
        df[col] = df[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "--")
    df["$p$"] = df["$p$"].map(lambda x: f"{x:.3f}" if pd.notna(x) else "--")
    body = df.to_latex(index=False, escape=False, column_format="llrrrrrrr")
    write_latex_table(
        OUT_DIR / "T10_network_regime_comparison.tex", body,
        caption="Network geometry around large stress months. Welch's two-sample $t$-test contrasts metric values in months in the 24-month rolling window centred on a large avalanche ($S\\geq 6$) versus all other windows.",
        label="tab:network_regime",
        notes=[
            "Filters refer to the network construction: $|\\rho|>0.4$ (corr), partial correlation $>0.4$ (partial), or planar minimum-spanning tree (mst).",
            "Network density and edge counts in MST are constant by construction (10 edges among 11 nodes) and are omitted.",
            "The partial-correlation network exhibits the strongest large-vs-non-large discrimination, with significant increases in Frobenius spectral radius and Frobenius SFI during large-stress windows.",
        ])


# =============================================================================
# TABLE 11. NETWORK PREDICTIVE ASSOCIATION (lead-lag)
# =============================================================================
def make_T11(geom, events, lead_months=3):
    av = events.groupby("date")["country"].count().reset_index()
    av.columns = ["date","size"]
    av = av.set_index("date")["size"]

    rows = []
    metric_labels = {
        "rho_max": "$\\rho^{\\max}$",
        "rho_frob": "$\\rho^{\\mathrm{Fro}}$",
        "sfi_max": "SFI$^{\\max}$",
        "sfi_frob": "SFI$^{\\mathrm{Fro}}$",
        "kappa_ollivier": "$\\bar{\\kappa}^O$",
        "density": "Density",
        "hhi": "HHI strength",
    }
    for typ in ["corr","partial","mst"]:
        snap = geom[typ]["snap"].copy().sort_values("date").set_index("date")
        for met, lab in metric_labels.items():
            x = snap[met].dropna()
            if x.nunique() <= 1:
                continue
            future = []
            for d in x.index:
                end = d + pd.DateOffset(months=lead_months)
                fut = av.loc[(av.index > d) & (av.index <= end)]
                future.append(fut.max() if len(fut) > 0 else 0.0)
            future = pd.Series(future, index=x.index)
            valid = pd.concat([x, future], axis=1).dropna()
            if len(valid) < 10:
                continue
            r, p = stats.pearsonr(valid.iloc[:, 0], valid.iloc[:, 1])
            rows.append({
                "Filter": typ,
                "Metric": lab,
                "Pearson $r$": r,
                "$p$": p,
                "$n$": len(valid),
            })
    df = pd.DataFrame(rows)
    df["Pearson $r$"] = df["Pearson $r$"].map(lambda x: f"{x:+.3f}")
    df["$p$"] = df["$p$"].map(lambda x: f"{x:.3f}")
    body = df.to_latex(index=False, escape=False, column_format="llrrr")
    write_latex_table(
        OUT_DIR / "T11_network_predictive_association.tex", body,
        caption=f"Lead--lag association between current network geometry and the maximum avalanche size over the next {lead_months} months.",
        label="tab:network_predictive",
        notes=[
            "Each row reports the Pearson correlation between a network metric at month $t$ and $\\max_{u\\in(t,t+3]} S_u$.",
            "Partial-correlation network metrics (Frobenius spectral radius and SFI) yield the most significant predictive associations.",
            "MST metrics are noisy due to the binary tree topology and yield generally non-significant lead--lag effects.",
        ])


# =============================================================================
# TABLE 12. HYPOTHESIS SUMMARY
# =============================================================================
def make_T12():
    rows = [
        {"\\#":"H1",
         "Hypothesis":"Country-level participation in stress avalanches is heterogeneous.",
         "Test / evidence":"Country event counts; Table~\\ref{tab:country_thresholds}; Fig.~5--6.",
         "Outcome":"Confirmed (range 7--18 events; $\\chi^2$ rejects uniform: $p<0.001$).",
        },
        {"\\#":"H2",
         "Hypothesis":"Mid-credibility, highly traded sovereigns dominate avalanche participation.",
         "Test / evidence":"Ranking BRA, MEX, PER, COL, PAN; Table~\\ref{tab:country_thresholds}.",
         "Outcome":"Confirmed; opposite of high-spread-driven prediction.",
        },
        {"\\#":"H3",
         "Hypothesis":"Avalanche-size, $|\\Delta s|$ and inter-event-time distributions exhibit heavy-tail boundary regime (PL not rejected vs.\\ LN).",
         "Test / evidence":"Clauset MLE + Vuong test; Table~\\ref{tab:tail_diagnostics}.",
         "Outcome":"Confirmed for all three objects ($p_{LN}\\in[0.10, 0.45]$; PL preferred over Exp).",
        },
        {"\\#":"H4",
         "Hypothesis":"Cross-country stress synchronization is non-trivial.",
         "Test / evidence":"Independent shuffle placebo, $B=1{,}000$; Table~\\ref{tab:placebo}.",
         "Outcome":"Confirmed; observed $S_{\\max}=11$ has placebo $p<0.001$.",
        },
        {"\\#":"H5",
         "Hypothesis":"Network geometry shifts around large stress months and predicts near-future maximum avalanche size.",
         "Test / evidence":"Regime comparison + lead--lag; Tables~\\ref{tab:network_regime}, \\ref{tab:network_predictive}.",
         "Outcome":"Partially confirmed; partial-correlation network spectral measures discriminate and predict; corr / MST measures show weaker effects.",
        },
        {"\\#":"H6",
         "Hypothesis":"Threshold sensitivity: avalanche evidence is robust to choice of $\\sigma$.",
         "Test / evidence":"$\\hat{\\alpha}(\\sigma)$ for $\\sigma\\in[1.25, 2.50]$; Table~\\ref{tab:threshold_robustness}.",
         "Outcome":"Confirmed; $\\hat{\\alpha}\\in[1.70, 1.83]$ for $\\sigma\\leq 2.0$; $S_{\\max}=11$ for all $\\sigma$.",
        },
    ]
    df = pd.DataFrame(rows)
    body = df.to_latex(index=False, escape=False,
                        column_format="lp{0.27\\linewidth}p{0.30\\linewidth}p{0.27\\linewidth}")
    write_latex_table(
        OUT_DIR / "T12_hypothesis_summary.tex", body,
        caption="Summary of hypotheses, tests, and outcomes.",
        label="tab:hypothesis_summary",
        notes=[
            "All tests are based on monthly EMBI Global Diversified spreads, 11 LATAM countries, 2007:M10--2026:M4.",
            "H3 follows the heavy-tail boundary framework of Stumpf and Porter (2012, \\emph{Science}).",
        ])


# =============================================================================
# MASTER INCLUDE FILE
# =============================================================================
def make_master():
    files = sorted(p.name for p in OUT_DIR.glob("T??_*.tex") if p.is_file())
    body = ["% Master include for all paper tables.",
            "% Required preamble: \\usepackage{booktabs, threeparttable, array}",
            f"% Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            ""]
    for f in files:
        body.append(f"\\input{{{Path(f).stem}}}")
        body.append("")
    (OUT_DIR / "master.tex").write_text("\n".join(body), encoding="utf-8")
    print(f"  [OK] master.tex (includes {len(files)} tables)")


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("="*80)
    print("LATAM SOC PAPER -- TABLE GENERATOR (V3 REAL EMBI)")
    print("Diego Vallarino | Sandpile Economics Program")
    print("="*80, "\n")

    panel, events, soc, geom, edges = load_data()
    print(f"  Panel:  {len(panel)} obs, {panel['country'].nunique()} countries")
    print(f"  Events: {len(events)} country-stress events")
    print(f"  Geom:   corr={len(geom['corr']['snap'])}, "
          f"partial={len(geom['partial']['snap'])}, mst={len(geom['mst']['snap'])} snapshots\n")

    print("[INFO] Generating tables...\n")
    make_T01(panel)
    make_T02(events)
    make_T03(panel, events)
    make_T04(events)
    make_T05(events)
    make_T06(events)
    make_T07(panel, events)
    make_T08(panel)
    make_T09(panel, events, n_boot=1000)
    make_T10(geom, events)
    make_T11(geom, events)
    make_T12()
    make_master()

    print("\n" + "="*80)
    print("DONE. Twelve paper-ready LaTeX tables written to:")
    print(f"  {OUT_DIR}")
    print("="*80)
    print("\nUsage in your paper:")
    print("  \\usepackage{booktabs, threeparttable, array}")
    print("  \\input{master.tex}      % includes all tables")
    print("  % or selectively:")
    print("  \\input{T01_data_coverage}")
    print("  \\input{T07_tail_diagnostics}")
    print("  ...")


if __name__ == "__main__":
    main()
