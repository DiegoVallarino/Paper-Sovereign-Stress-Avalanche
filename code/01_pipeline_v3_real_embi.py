"""
================================================================================
LATAM FISCAL VULNERABILITY: SOC AND NETWORK GEOMETRY
                          VERSION 3.0 — REAL EMBI DATA
================================================================================
Author: Diego Vallarino, Ph.D.
Affiliation: Inter-American Development Bank
Project: Sandpile Economics Program

CHANGES vs. V2.1:
  - Replaces the synthetic sandpile generator with REAL JPM EMBI Global
    Diversified spreads, monthly, 2007-10 to 2026-04.
  - 11 LATAM countries: ARG, BRA, CHL, COL, ECU, MEX, PAN, PER, DOM, URY, SLV.
  - Empirical avalanches built from spread jumps exceeding country-specific
    1.5-sigma thresholds on dlog(spread).
  - Sandpile model is REPOSITIONED as a structural-mechanism appendix, not as
    primary data source. The empirical paper uses real spreads end-to-end.
  - Network geometry (Forman + Ollivier curvature, spectral radius) computed
    on rolling 24-month correlation windows of REAL spread changes.
  - Same statistical framework as V2.1 (Clauset xmin, lognormal comparison).

INPUT:
  - /mnt/user-data/uploads/EMBI_mensual.xlsx  (or local copy)
  - World Bank fiscal indicators (used only for country covariates in
    secondary analyses, NOT to drive the dynamics).

OUTPUT:
  - /home/diego/Escritorio/figures/ (or Windows path on Escritorio)
  - /home/diego/Escritorio/data_latam_fiscal/

PAPER NARRATIVE:
  This is now an EMPIRICAL paper. Title candidates:
    1. "Endogenous Fragility in LATAM Sovereign Markets: Evidence from
        EMBI Avalanches and Network Geometry, 2007-2026"
    2. "Self-Organized Criticality in Latin American Sovereign Risk:
        Empirical Evidence from JPM EMBI Spreads"
    3. "The Geometry of Sovereign Contagion: A Network-Curvature Analysis
        of LATAM EMBI Spreads"
================================================================================
"""

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats
from scipy.linalg import eigvals
from scipy.optimize import linprog

warnings.filterwarnings("ignore")

# ==============================================================================
# CONFIG
# ==============================================================================

if sys.platform.startswith("win"):
    BASE_DIR = Path(r"C:\Users\diego\OneDrive\Escritorio\figures")
    EMBI_PATH = Path(r"C:\Users\diego\OneDrive\Escritorio\EMBI_mensual.xlsx")
else:
    BASE_DIR = Path.home() / "Escritorio" / "figures"
    EMBI_PATH = Path.home() / "Escritorio" / "EMBI_mensual.xlsx"

FIG_DIR = BASE_DIR
DATA_DIR = BASE_DIR.parent / "data_latam_fiscal"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

print(f"[INFO] EMBI source: {EMBI_PATH}")
print(f"[INFO] Figures:     {FIG_DIR}")
print(f"[INFO] Data:        {DATA_DIR}")

LATAM_COUNTRIES = {
    "ARG":"Argentina","BRA":"Brazil","CHL":"Chile","COL":"Colombia",
    "ECU":"Ecuador","MEX":"Mexico","PAN":"Panama","PER":"Peru",
    "DOM":"Dominican Rep.","URY":"Uruguay","SLV":"El Salvador",
}
# Stress threshold in standard deviations of dlog(spread)
STRESS_SIGMA = 1.5

sns.set_style("whitegrid")
plt.rcParams.update({
    "font.size":11,"axes.labelsize":11,"axes.titlesize":12,
    "legend.fontsize":9,"figure.dpi":110,"savefig.dpi":300,
    "savefig.bbox":"tight","font.family":"serif",
})


# ==============================================================================
# 1. LOAD REAL EMBI DATA (replaces sandpile generator)
# ==============================================================================

def load_embi_panel(filepath=EMBI_PATH):
    """
    Load JPM EMBI Global Diversified monthly stripped spreads.
    Input file: Excel with daily and monthly blocks; we use the monthly block.
    Spreads come in percentage points; we convert to bps (multiply by 100).
    Returns long panel: date, country, spread_bp.
    """
    print(f"[INFO] Loading EMBI panel from {filepath}...")
    df = pd.read_excel(filepath, sheet_name="Hoja1", header=None)

    # Header at row 7, monthly block at columns 19-32
    data = df.iloc[8:].copy()
    monthly = data[[19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]].copy()
    monthly.columns = ["date","month","year","ARG","BRA","CHL","COL","ECU",
                       "MEX","PAN","PER","DOM","URY","SLV"]
    monthly["date"] = pd.to_datetime(monthly["date"], errors="coerce")
    monthly = monthly.dropna(subset=["date"]).reset_index(drop=True)

    countries = list(LATAM_COUNTRIES.keys())
    for c in countries:
        monthly[c] = pd.to_numeric(monthly[c], errors="coerce") * 100  # to bps

    # Long format
    rows = []
    for _, r in monthly.iterrows():
        for c in countries:
            if pd.notna(r[c]):
                rows.append({"date": r["date"], "country": c,
                             "spread_bp": float(r[c])})
    panel = pd.DataFrame(rows).sort_values(["country", "date"]).reset_index(drop=True)
    panel.to_csv(DATA_DIR / "latam_spread_panel_real_embi.csv", index=False)
    print(f"  [OK] {len(panel)} obs, {panel['country'].nunique()} countries, "
          f"{panel['date'].min().strftime('%Y-%m')} to {panel['date'].max().strftime('%Y-%m')}")
    return panel


def build_empirical_avalanches(spread_panel, sigma=STRESS_SIGMA):
    """
    Empirical avalanche definition:
    A country is hit in month t if its dlog(spread) exceeds sigma*sd of its
    own dlog series. An avalanche of size S is a month with S countries hit.
    """
    print(f"[INFO] Building empirical avalanches (sigma={sigma})...")
    rows = []
    countries = sorted(spread_panel["country"].unique())
    for c in countries:
        s = spread_panel[spread_panel["country"]==c].sort_values("date").set_index("date")
        dlog = np.log(s["spread_bp"]).diff()
        sd = dlog.std()
        if pd.isna(sd) or sd == 0:
            continue
        thresh = sigma * sd
        for d, v in dlog.items():
            if pd.notna(v) and v > thresh:
                rows.append({
                    "date": d, "country": c,
                    "dlog_spread": float(v),
                    "overshoot": float(v - thresh),
                    "jump_bp": float(s["spread_bp"].diff().loc[d]),
                    "spread_bp": float(s["spread_bp"].loc[d]),
                })
    event_df = pd.DataFrame(rows).sort_values("date").reset_index(drop=True)
    event_df.to_csv(DATA_DIR / "latam_avalanche_events_real.csv", index=False)
    print(f"  [OK] {len(event_df)} country-stress events identified")
    if len(event_df) > 0:
        top = event_df["country"].value_counts().head()
        print(f"  Top stressed: " + ", ".join([f"{c}({n})" for c,n in top.items()]))
    return event_df


# ==============================================================================
# 2. POWER-LAW ESTIMATION (same as V2.1)
# ==============================================================================

def estimate_powerlaw_clauset(data, min_n=20, discrete=False):
    try:
        import powerlaw
    except ImportError:
        return {"alpha":np.nan,"xmin":np.nan,"ks_D":np.nan,"p_pl":np.nan,
                "R_lognormal":np.nan,"p_lognormal":np.nan,"n":0,"n_tail":0}
    x = np.asarray(data, dtype=float)
    x = x[(x>0) & np.isfinite(x)]
    if len(x) < min_n:
        return {"alpha":np.nan,"xmin":np.nan,"ks_D":np.nan,"p_pl":np.nan,
                "R_lognormal":np.nan,"p_lognormal":np.nan,
                "n":int(len(x)),"n_tail":0}
    try:
        fit = powerlaw.Fit(x, discrete=discrete, verbose=False)
        alpha = float(fit.power_law.alpha)
        xmin = float(fit.power_law.xmin)
        D = float(fit.power_law.D)
        try:
            R, p_ln = fit.distribution_compare("power_law","lognormal",
                                               normalized_ratio=True)
        except Exception:
            R, p_ln = np.nan, np.nan
        n_tail = int((x >= xmin).sum())
        return {"alpha":alpha,"xmin":xmin,"ks_D":D,"p_pl":np.nan,
                "R_lognormal":float(R) if pd.notna(R) else np.nan,
                "p_lognormal":float(p_ln) if pd.notna(p_ln) else np.nan,
                "n":int(len(x)),"n_tail":n_tail}
    except Exception:
        return {"alpha":np.nan,"xmin":np.nan,"ks_D":np.nan,"p_pl":np.nan,
                "R_lognormal":np.nan,"p_lognormal":np.nan,
                "n":int(len(x)),"n_tail":0}


# ==============================================================================
# 3. INTER-EVENT TIMES + COUNTRY-LEVEL SOC
# ==============================================================================

def compute_country_iet_soc(event_df, dates, window_months=60):
    print(f"[INFO] Country-level IET SOC (window={window_months}m)...")
    if len(event_df) == 0:
        return pd.DataFrame()
    rows = []
    for c in sorted(event_df["country"].unique()):
        ev_c = event_df[event_df["country"]==c].sort_values("date").reset_index(drop=True)
        if len(ev_c) < 3:
            continue
        ev_c["date"] = pd.to_datetime(ev_c["date"])
        for end_date in pd.to_datetime(dates):
            start_date = end_date - pd.DateOffset(months=window_months)
            mask = (ev_c["date"] >= start_date) & (ev_c["date"] <= end_date)
            ev_w = ev_c[mask].sort_values("date").reset_index(drop=True)
            if len(ev_w) < 3:
                continue
            wait = ev_w["date"].diff().dt.days.dropna().values / 30.44
            wait = wait[wait > 0]
            if len(wait) < 3:
                continue
            res = estimate_powerlaw_clauset(wait, min_n=3)
            rows.append({
                "date": end_date, "country": c,
                "alpha_iet": res["alpha"], "xmin_iet": res["xmin"],
                "n_events_window": len(ev_w),
                "mean_wait_months": float(np.mean(wait)),
            })
    df = pd.DataFrame(rows)
    df.to_csv(DATA_DIR / "soc_intereventimes_panel_real.csv", index=False)
    valid = df["alpha_iet"].notna().sum() if len(df) > 0 else 0
    print(f"  [OK] {len(df)} obs (alpha valid: {valid})")
    return df


# ==============================================================================
# 4. CONTAGION NETWORKS (same logic, on REAL spreads)
# ==============================================================================

def partial_correlation_matrix(X):
    common = X.mean(axis=1)
    resid = X.sub(common, axis=0)
    return resid.corr()


def build_contagion_networks(spread_panel, window_months=24, step_months=3,
                              corr_threshold=0.40):
    print(f"[INFO] Networks (window={window_months}m, corr>{corr_threshold})...")
    pivot = spread_panel.pivot(index="date", columns="country",
                                values="spread_bp").sort_index()
    dspread = pivot.diff()
    countries = list(pivot.columns)
    nets_corr, nets_partial, nets_mst = {}, {}, {}
    edges = []
    dates = pivot.index
    for end_idx in range(window_months, len(dates), step_months):
        start_idx = end_idx - window_months
        wnd = dspread.iloc[start_idx:end_idx]
        if wnd.dropna(how="all").shape[0] < window_months // 2:
            continue
        d = dates[end_idx]

        corr = wnd.corr().abs()
        Gc = nx.Graph(); [Gc.add_node(c) for c in countries]
        for i,c1 in enumerate(countries):
            for c2 in countries[i+1:]:
                w = corr.loc[c1,c2]
                if pd.notna(w) and w > corr_threshold:
                    Gc.add_edge(c1,c2,weight=float(w))
                    edges.append({"date":d,"type":"corr","src":c1,"dst":c2,"weight":float(w)})
        nets_corr[d] = Gc

        try:
            pcorr = partial_correlation_matrix(wnd.dropna()).abs()
        except Exception:
            pcorr = corr
        Gp = nx.Graph(); [Gp.add_node(c) for c in countries]
        for i,c1 in enumerate(countries):
            for c2 in countries[i+1:]:
                w = pcorr.loc[c1,c2]
                if pd.notna(w) and w > corr_threshold:
                    Gp.add_edge(c1,c2,weight=float(w))
                    edges.append({"date":d,"type":"partial","src":c1,"dst":c2,"weight":float(w)})
        nets_partial[d] = Gp

        Gfull = nx.Graph(); [Gfull.add_node(c) for c in countries]
        for i,c1 in enumerate(countries):
            for c2 in countries[i+1:]:
                w = corr.loc[c1,c2]
                if pd.notna(w) and w > 0:
                    Gfull.add_edge(c1,c2,distance=float(1-w),weight=float(w))
        try:
            T_mst = nx.minimum_spanning_tree(Gfull, weight="distance")
        except Exception:
            T_mst = nx.Graph(); [T_mst.add_node(c) for c in countries]
        nets_mst[d] = T_mst
        for u,v,a in T_mst.edges(data=True):
            edges.append({"date":d,"type":"mst","src":u,"dst":v,"weight":float(a.get("weight",0))})

    pd.DataFrame(edges).to_csv(DATA_DIR/"fiscal_network_edges_real.csv", index=False)
    print(f"  [OK] corr={len(nets_corr)}, partial={len(nets_partial)}, mst={len(nets_mst)}")
    return nets_corr, nets_partial, nets_mst


# ==============================================================================
# 5. NETWORK GEOMETRY (with FIXED Ollivier-Ricci)
# ==============================================================================

def spectral_metrics(G, nodelist=None):
    if nodelist is None:
        nodelist = list(G.nodes())
    if G.number_of_edges() == 0:
        return {"rho_max":0.0,"rho_frob":0.0,"sfi_max":0.0,"sfi_frob":0.0}
    A = nx.to_numpy_array(G, nodelist=nodelist, weight="weight")
    rs = A.sum(axis=1)
    max_rs = float(rs.max()) if rs.size > 0 else 0.0
    fro = float(np.sqrt(np.sum(A*A)))
    eigs = np.abs(eigvals(A))
    lambda_max = float(np.max(eigs)) if eigs.size > 0 else 0.0
    rho_max = lambda_max / max_rs if max_rs > 1e-10 else 0.0
    rho_frob = lambda_max / fro if fro > 1e-10 else 0.0
    rho_max = min(rho_max, 0.9999)
    rho_frob = min(rho_frob, 0.9999)
    sfi_max = rho_max / (1 - rho_max) if rho_max < 1 else np.nan
    sfi_frob = rho_frob / (1 - rho_frob) if rho_frob < 1 else np.nan
    return {"rho_max":rho_max,"rho_frob":rho_frob,
            "sfi_max":sfi_max,"sfi_frob":sfi_frob}


def ollivier_ricci_curvature(G, edge):
    """
    FIXED: uses unit baseline distance (Ollivier original), bounded in [-1, 1].
    """
    u, v = edge
    if not G.has_edge(u, v):
        return np.nan
    Nu = list(G.neighbors(u))
    Nv = list(G.neighbors(v))
    if not Nu or not Nv:
        return 0.0
    try:
        # Distance graph: 1 - normalized weight. Use shortest path.
        max_w = max((d.get("weight",1) for _,_,d in G.edges(data=True)), default=1.0)
        D = nx.Graph()
        for nn in G.nodes():
            D.add_node(nn)
        for x, y, dat in G.edges(data=True):
            w = dat.get("weight", 1.0)
            D.add_edge(x, y, distance=max(0.01, 1.0 - w/max_w))

        mu_u = {n: 1.0/len(Nu) for n in Nu}
        mu_v = {n: 1.0/len(Nv) for n in Nv}

        dist_uv = {}
        for nu in Nu:
            for nv in Nv:
                try:
                    dist_uv[(nu, nv)] = nx.shortest_path_length(D, nu, nv, weight="distance")
                except nx.NetworkXNoPath:
                    dist_uv[(nu, nv)] = 5.0

        nu_n, nv_n = len(Nu), len(Nv)
        c = np.array([dist_uv[(Nu[i], Nv[j])]
                      for i in range(nu_n) for j in range(nv_n)])
        A_eq = np.zeros((nu_n + nv_n, nu_n*nv_n))
        b_eq = np.zeros(nu_n + nv_n)
        for i in range(nu_n):
            for j in range(nv_n):
                A_eq[i, i*nv_n + j] = 1.0
            b_eq[i] = mu_u[Nu[i]]
        for j in range(nv_n):
            for i in range(nu_n):
                A_eq[nu_n + j, i*nv_n + j] = 1.0
            b_eq[nu_n + j] = mu_v[Nv[j]]
        bounds = [(0, None)] * (nu_n*nv_n)
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
        if not res.success:
            return np.nan
        W1 = float(res.fun)
        # FIX: use unit baseline (Ollivier original). Bounded in [-1, 1] approximately.
        kappa = 1.0 - W1
        return float(kappa)
    except Exception:
        return np.nan


def forman_ricci_curvature(G, edge):
    u, v = edge
    if not G.has_edge(u, v):
        return np.nan
    we = G[u][v]["weight"]
    if we <= 0:
        return np.nan
    wu = sum(G[u][k]["weight"] for k in G.neighbors(u))
    wv = sum(G[v][k]["weight"] for k in G.neighbors(v))
    s_u = sum(we / np.sqrt(we * G[u][k]["weight"])
              for k in G.neighbors(u) if k != v and G[u][k]["weight"] > 0)
    s_v = sum(we / np.sqrt(we * G[v][k]["weight"])
              for k in G.neighbors(v) if k != u and G[v][k]["weight"] > 0)
    return float(we * ((wu + wv) / we - s_u - s_v))


def compute_network_geometry(networks, label="corr"):
    print(f"[INFO] Geometry for '{label}'...")
    snap_rows, country_rows = [], []
    for date, G in networks.items():
        if G.number_of_edges() == 0:
            continue
        kappa_O = [ollivier_ricci_curvature(G, e) for e in G.edges()]
        kappa_O = [c for c in kappa_O if pd.notna(c)]
        kappa_F = [forman_ricci_curvature(G, e) for e in G.edges()]
        kappa_F = [c for c in kappa_F if pd.notna(c)]
        kappa_O_bar = float(np.mean(kappa_O)) if kappa_O else np.nan
        kappa_F_bar = float(np.mean(kappa_F)) if kappa_F else np.nan
        spec = spectral_metrics(G)
        nodes = list(G.nodes())
        strengths = np.array([sum(G[u][k]["weight"] for k in G.neighbors(u))
                              for u in nodes])
        if strengths.sum() > 0:
            shares = strengths / strengths.sum()
            hhi = float(np.sum(shares**2))
        else:
            hhi = np.nan
        density = float(nx.density(G))
        snap_rows.append({
            "date":date,"type":label,
            "kappa_ollivier":kappa_O_bar,"kappa_forman":kappa_F_bar,
            "rho_max":spec["rho_max"],"rho_frob":spec["rho_frob"],
            "sfi_max":spec["sfi_max"],"sfi_frob":spec["sfi_frob"],
            "hhi":hhi,"density":density,
            "n_nodes":G.number_of_nodes(),"n_edges":G.number_of_edges(),
        })
        for u in nodes:
            edges_u = list(G.edges(u))
            if not edges_u:
                country_rows.append({"date":date,"type":label,"country":u,
                                     "kappa_country":np.nan,
                                     "strength":0.0,"degree":0})
                continue
            curvs_u = [ollivier_ricci_curvature(G, e) for e in edges_u]
            curvs_u = [c for c in curvs_u if pd.notna(c)]
            country_rows.append({
                "date":date,"type":label,"country":u,
                "kappa_country":float(np.mean(curvs_u)) if curvs_u else np.nan,
                "strength":float(strengths[nodes.index(u)]),
                "degree":int(G.degree(u)),
            })
    snap_df = pd.DataFrame(snap_rows).sort_values("date").reset_index(drop=True)
    country_df = pd.DataFrame(country_rows).sort_values(["country","date"]).reset_index(drop=True)
    snap_df.to_csv(DATA_DIR/f"network_geometry_snapshots_{label}_real.csv", index=False)
    country_df.to_csv(DATA_DIR/f"network_geometry_country_{label}_real.csv", index=False)
    print(f"  [OK] {len(snap_df)} snapshots, {len(country_df)} country-snap obs")
    return snap_df, country_df


# ==============================================================================
# 6. FIGURES
# ==============================================================================

def fig1_real_spreads(spread_panel, event_df):
    """Figure 1: real EMBI spreads with avalanche dates highlighted."""
    fig, ax = plt.subplots(figsize=(12, 6))
    pivot = spread_panel.pivot(index="date", columns="country", values="spread_bp")
    palette = sns.color_palette("husl", len(pivot.columns))
    for c, col in zip(pivot.columns, palette):
        ax.plot(pivot.index, pivot[c], lw=1.0, alpha=0.75, label=c, color=col)

    # Highlight big avalanches
    if len(event_df) > 0:
        av = event_df.groupby("date")["country"].count().reset_index()
        av.columns = ["date", "size"]
        big = av[av["size"] >= 6]
        for _, r in big.iterrows():
            ax.axvline(r["date"], color="red", alpha=0.25, lw=0.8)

    ax.set_yscale("log")
    ax.set_ylabel("EMBI spread (bps, log)")
    ax.set_xlabel("Year")
    ax.set_title("Figure 1. LATAM EMBI spreads (real data, 2007–2026); "
                 "red lines = months with $\\geq$6 countries in stress")
    ax.legend(ncol=6, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.10))
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig1_latam_spreads.png")
    plt.savefig(FIG_DIR/"fig1_latam_spreads.pdf")
    plt.close()
    print(f"  [OK] fig1")


def fig2_avalanche_powerlaw(event_df):
    if len(event_df) == 0:
        return
    av = event_df.groupby("date")["country"].count().reset_index()
    av.columns = ["date", "size"]
    sizes = av["size"].values

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: avalanche CCDF
    res = estimate_powerlaw_clauset(sizes, min_n=10, discrete=True)
    s_sorted = np.sort(sizes)[::-1]
    ranks = np.arange(1, len(s_sorted)+1) / len(s_sorted)
    axes[0].loglog(s_sorted, ranks, 'o', ms=7, color="darkred",
                   alpha=0.7, label="Empirical")
    if pd.notna(res["alpha"]) and pd.notna(res["xmin"]):
        x_fit = np.linspace(res["xmin"], s_sorted.max(), 80)
        ccdf_fit = (x_fit / res["xmin"]) ** -(res["alpha"] - 1)
        ccdf_fit *= (sizes >= res["xmin"]).mean()
        axes[0].loglog(x_fit, ccdf_fit, "k--", lw=2,
                       label=f"PL fit α̂={res['alpha']:.2f}")
    axes[0].set_xlabel("Avalanche size (countries hit per month)")
    axes[0].set_ylabel("CCDF P(S > s)")
    axes[0].set_title("Panel A. Avalanche size CCDF (real EMBI)")
    info = (f"α̂ = {res['alpha']:.2f}\n"
            f"x_min = {res['xmin']:.1f}\n"
            f"R vs lognorm = {res['R_lognormal']:+.2f}\n"
            f"p_LN = {res['p_lognormal']:.3f}\n"
            f"n_events = {len(sizes)}")
    axes[0].text(0.55, 0.65, info, transform=axes[0].transAxes,
                 fontsize=9, va="top",
                 bbox=dict(boxstyle="round", fc="white", alpha=0.9))
    axes[0].legend(loc="lower left", fontsize=9)
    axes[0].grid(True, which="both", alpha=0.3)

    # Panel B: events per country
    ec = event_df["country"].value_counts().sort_values()
    palette = sns.color_palette("RdYlGn_r", len(ec))
    axes[1].barh(ec.index, ec.values, color=palette, edgecolor="black",
                 linewidth=0.4)
    axes[1].set_xlabel("Number of stress events")
    axes[1].set_title("Panel B. Country participation in stress events")
    axes[1].grid(True, alpha=0.3, axis="x")

    plt.suptitle("Figure 2. Empirical avalanche dynamics in LATAM EMBI spreads",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig2_avalanche_soc.png")
    plt.savefig(FIG_DIR/"fig2_avalanche_soc.pdf")
    plt.close()
    print(f"  [OK] fig2")


def fig3_networks(nets_corr, nets_partial, nets_mst, dates_to_plot=None):
    if dates_to_plot is None:
        all_dates = sorted(nets_corr.keys())
        # Pick four meaningful dates
        target_dates = [pd.Timestamp("2009-01-01"),  # GFC
                        pd.Timestamp("2013-01-01"),  # Taper tantrum aftermath
                        pd.Timestamp("2018-04-01"),  # Argentina crisis
                        pd.Timestamp("2021-04-01")]  # Post-COVID
        dates_to_plot = []
        for td in target_dates:
            closest = min(all_dates, key=lambda d: abs(pd.Timestamp(d) - td))
            dates_to_plot.append(closest)

    fig, axes = plt.subplots(3, len(dates_to_plot),
                             figsize=(5*len(dates_to_plot), 14))
    if len(dates_to_plot) == 1:
        axes = axes.reshape(3, 1)
    sample_G = list(nets_corr.values())[-1]
    pos = nx.spring_layout(sample_G, seed=42, k=0.7)

    for col, d in enumerate(dates_to_plot):
        for row, (label, nets) in enumerate([("Corr",nets_corr),
                                              ("Partial",nets_partial),
                                              ("MST",nets_mst)]):
            ax = axes[row, col]
            G = nets.get(d, nx.Graph())
            if G.number_of_edges() == 0:
                ax.set_title(f"{label}: {pd.Timestamp(d).strftime('%Y-%m')} (empty)")
                ax.axis("off")
                continue
            weights = [G[u][v].get("weight",0)*4 for u,v in G.edges()]
            strengths = [sum(G[u][k].get("weight",0) for k in G.neighbors(u))
                         for u in G.nodes()]
            sizes = [200 + 100*s for s in strengths]
            nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes,
                                    node_color="#1f77b4", alpha=0.85,
                                    edgecolors="black", linewidths=0.6)
            nx.draw_networkx_edges(G, pos, ax=ax, width=weights,
                                    alpha=0.55, edge_color="gray")
            nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
            ax.set_title(f"{label}: {pd.Timestamp(d).strftime('%Y-%m')}\n"
                         f"edges={G.number_of_edges()}, "
                         f"density={nx.density(G):.2f}", fontsize=9)
            ax.axis("off")

    plt.suptitle("Figure 3. LATAM contagion networks (real EMBI): three filtering methods",
                 fontsize=13, y=1.00)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig3_network_snapshots.png")
    plt.savefig(FIG_DIR/"fig3_network_snapshots.pdf")
    plt.close()
    print(f"  [OK] fig3")


def fig4_geometry(snap_dfs):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    palette = {"corr":"darkred","partial":"navy","mst":"darkgreen"}
    for snap_df, label in snap_dfs:
        if len(snap_df) == 0:
            continue
        col = palette.get(label, "black")
        axes[0,0].plot(snap_df["date"], snap_df["kappa_ollivier"], color=col, lw=1.3, label=label)
        axes[0,1].plot(snap_df["date"], snap_df["kappa_forman"], color=col, lw=1.3, label=label)
        axes[1,0].plot(snap_df["date"], snap_df["rho_max"], color=col, lw=1.3, label=label)
        axes[1,1].plot(snap_df["date"], snap_df["rho_frob"], color=col, lw=1.3, label=label)
        axes[2,0].plot(snap_df["date"], snap_df["sfi_max"], color=col, lw=1.3, label=label)
        axes[2,1].plot(snap_df["date"], snap_df["density"], color=col, lw=1.3, label=label)
    axes[0,0].set_ylabel(r"$\bar{\kappa}^{O}_t$"); axes[0,0].set_title("Panel A. Ollivier-Ricci curvature (main)")
    axes[0,1].set_ylabel(r"$\bar{\kappa}^{F}_t$"); axes[0,1].set_title("Panel B. Forman-Ricci curvature (auxiliary)")
    axes[1,0].set_ylabel(r"$\rho^{\max}_t$"); axes[1,0].set_title("Panel C. Spectral radius (max-row-sum)")
    axes[1,1].set_ylabel(r"$\rho^{\mathrm{Fro}}_t$"); axes[1,1].set_title("Panel D. Spectral radius (Frobenius)")
    axes[2,0].set_ylabel(r"SFI"); axes[2,0].set_title("Panel E. Spectral fragility index")
    axes[2,1].set_ylabel("Density"); axes[2,1].set_title("Panel F. Network density")
    for ax in axes.flat:
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    for ax in axes[2, :]:
        ax.set_xlabel("Year")
    plt.suptitle("Figure 4. Network geometry on real EMBI (V3.0)", fontsize=13, y=1.00)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig4_geometry_dynamics.png")
    plt.savefig(FIG_DIR/"fig4_geometry_dynamics.pdf")
    plt.close()
    print(f"  [OK] fig4")


def fig5_alpha_iet(soc_iet):
    if len(soc_iet) == 0 or soc_iet["alpha_iet"].notna().sum() == 0:
        return
    valid = soc_iet.dropna(subset=["alpha_iet"]).copy()
    fig, ax = plt.subplots(figsize=(11, 5.8))
    palette = sns.color_palette("husl", valid["country"].nunique())
    for (c, g), col in zip(valid.groupby("country"), palette):
        ax.plot(g["date"], g["alpha_iet"], lw=1.2, alpha=0.8, label=c, color=col)
    ax.axhline(2.0, color="red", ls="--", lw=1, label=r"$\alpha=2$")
    ax.axhline(3.0, color="orange", ls="--", lw=1, label=r"$\alpha=3$")
    ax.set_ylabel(r"$\hat{\alpha}^{IET}_{c,t}$ (rolling 60m)")
    ax.set_xlabel("Year")
    ax.set_title("Figure 5. Country-level SOC: rolling tail exponent of inter-event times (real EMBI)")
    ax.legend(ncol=6, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.10))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig5_soc_country_panel.png")
    plt.savefig(FIG_DIR/"fig5_soc_country_panel.pdf")
    plt.close()
    print(f"  [OK] fig5")


def fig6_alpha_kappa(soc_iet, country_corr):
    if len(soc_iet) == 0 or len(country_corr) == 0:
        return
    soc_q = soc_iet.dropna(subset=["alpha_iet"]).copy()
    if len(soc_q) == 0:
        return
    soc_q["date"] = pd.to_datetime(soc_q["date"])
    soc_q["q"] = soc_q["date"].dt.to_period("Q").dt.to_timestamp()
    geom_q = country_corr.copy()
    geom_q["date"] = pd.to_datetime(geom_q["date"])
    geom_q["q"] = geom_q["date"].dt.to_period("Q").dt.to_timestamp()
    soc_qg = soc_q.groupby(["country","q"])[["alpha_iet"]].mean().reset_index()
    geom_qg = geom_q.groupby(["country","q"])[["kappa_country"]].mean().reset_index()
    merged = soc_qg.merge(geom_qg, on=["country","q"]).dropna()
    if len(merged) == 0:
        return
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sns.scatterplot(data=merged, x="kappa_country", y="alpha_iet",
                    hue="country", alpha=0.7, s=30, ax=ax, legend="full")
    ax.set_xlabel(r"Country mean Ollivier-Ricci curvature $\bar{\kappa}^{O}_{c,t}$")
    ax.set_ylabel(r"Tail exponent $\hat{\alpha}^{IET}_{c,t}$")
    ax.set_title("Figure 6. Joint behavior of country SOC and network geometry (real EMBI)")
    if len(merged) > 5:
        slope, intercept, r, p, _ = stats.linregress(merged["kappa_country"], merged["alpha_iet"])
        xs = np.linspace(merged["kappa_country"].min(), merged["kappa_country"].max(), 80)
        ax.plot(xs, slope*xs + intercept, "k--", lw=1.3,
                label=f"OLS r={r:.3f}, p={p:.3f}")
    ax.legend(fontsize=8, ncol=2, loc="best")
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig6_alpha_vs_kappa.png")
    plt.savefig(FIG_DIR/"fig6_alpha_vs_kappa.pdf")
    plt.close()
    print(f"  [OK] fig6")


def fig7_avalanche_timeline(event_df):
    if len(event_df) == 0:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    av = event_df.groupby("date")["country"].count().reset_index()
    av.columns = ["date", "size"]
    sizes = av["size"].values

    if len(sizes) > 5:
        s_sorted = np.sort(sizes)[::-1]
        ranks = np.arange(1, len(s_sorted)+1) / len(s_sorted)
        axes[0].loglog(s_sorted, ranks, 'o', ms=6, color="darkred", alpha=0.7)
        res = estimate_powerlaw_clauset(sizes, min_n=10, discrete=True)
        if pd.notna(res["alpha"]):
            axes[0].set_title(f"Panel A. Avalanche CCDF (α̂={res['alpha']:.2f})")
        else:
            axes[0].set_title("Panel A. Avalanche CCDF")
        axes[0].set_xlabel("Avalanche size")
        axes[0].set_ylabel("CCDF")
        axes[0].grid(True, which="both", alpha=0.3)

    ec = event_df.copy()
    ec["year"] = pd.to_datetime(ec["date"]).dt.year
    pivot = ec.pivot_table(index="year", columns="country",
                            values="overshoot", aggfunc="count", fill_value=0)
    pivot.plot(kind="bar", stacked=True, ax=axes[1],
               colormap="tab20", width=0.85, legend=False)
    axes[1].set_ylabel("Stress events per year")
    axes[1].set_xlabel("Year")
    axes[1].set_title("Panel B. Temporal distribution of stress events")
    axes[1].legend(ncol=6, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    plt.suptitle("Figure 7. Avalanche size and timing (real EMBI)", fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig7_avalanches.png")
    plt.savefig(FIG_DIR/"fig7_avalanches.pdf")
    plt.close()
    print(f"  [OK] fig7")


def fig8_dspread_diagnostic(spread_panel):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    all_d = []
    for c, g in spread_panel.groupby("country"):
        all_d.extend(g.sort_values("date")["spread_bp"].diff().dropna().values)
    x = np.abs(np.array(all_d))
    x = x[x > 0]
    res = estimate_powerlaw_clauset(x, min_n=50)
    x_sorted = np.sort(x)[::-1]
    ranks = np.arange(1, len(x_sorted)+1) / len(x_sorted)
    axes[0].loglog(x_sorted, ranks, 'o', ms=2, color="navy", alpha=0.5, label="Empirical")
    if pd.notna(res["alpha"]):
        x_fit = np.logspace(np.log10(res["xmin"]), np.log10(x_sorted.max()), 80)
        ccdf_pl = (x_fit / res["xmin"]) ** -(res["alpha"]-1)
        ccdf_pl *= (x >= res["xmin"]).mean()
        axes[0].loglog(x_fit, ccdf_pl, "r--", lw=1.5,
                       label=f"PL (α̂={res['alpha']:.2f})")
    axes[0].set_xlabel(r"$|\Delta s|$ (bps)")
    axes[0].set_ylabel("CCDF")
    axes[0].set_title("Panel A. |Δspread| pooled — heavy tail")
    info = (f"α̂ = {res['alpha']:.2f}\n"
            f"x_min = {res['xmin']:.1f}\n"
            f"R vs LN = {res['R_lognormal']:+.2f}\n"
            f"p_LN = {res['p_lognormal']:.3f}\n"
            f"(R<0 LN; |R|<2 ambiguous)")
    axes[0].text(0.04, 0.04, info, transform=axes[0].transAxes,
                 fontsize=9, va="bottom",
                 bbox=dict(boxstyle="round", fc="white", alpha=0.9))
    axes[0].legend(loc="upper right", fontsize=9)
    axes[0].grid(True, which="both", alpha=0.3)

    log_x = np.log(x[x > 1])
    axes[1].hist(log_x, bins=60, density=True, color="navy", alpha=0.6)
    mu, sd = log_x.mean(), log_x.std()
    xs = np.linspace(log_x.min(), log_x.max(), 200)
    axes[1].plot(xs, stats.norm.pdf(xs, mu, sd), "r-", lw=2,
                 label=f"N(μ={mu:.2f}, σ={sd:.2f})")
    axes[1].set_xlabel(r"$\log|\Delta s|$")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Panel B. log|Δs| histogram")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)
    plt.suptitle("Figure 8. Statistical diagnostic on real |Δspread| (EMBI)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR/"fig8_dspread_diagnostic.png")
    plt.savefig(FIG_DIR/"fig8_dspread_diagnostic.pdf")
    plt.close()
    print(f"  [OK] fig8")


# ==============================================================================
# 7. MAIN
# ==============================================================================

def main():
    print("=" * 80)
    print("LATAM SOC + NETWORK PIPELINE -- VERSION 3.0 (REAL EMBI DATA)")
    print("Diego Vallarino | Sandpile Economics Program")
    print("=" * 80)

    # 1. Load real EMBI
    spread_panel = load_embi_panel()

    # 2. Build empirical avalanches
    event_df = build_empirical_avalanches(spread_panel, sigma=STRESS_SIGMA)

    # 3. Country-level SOC from inter-event times
    dates = sorted(spread_panel["date"].unique())
    soc_iet = compute_country_iet_soc(event_df, dates, window_months=60)

    # 4. Contagion networks (correlation, partial, MST)
    nets_corr, nets_partial, nets_mst = build_contagion_networks(
        spread_panel, window_months=24, step_months=3, corr_threshold=0.40)

    # 5. Geometry
    snap_corr, country_corr = compute_network_geometry(nets_corr, "corr")
    snap_partial, country_partial = compute_network_geometry(nets_partial, "partial")
    snap_mst, country_mst = compute_network_geometry(nets_mst, "mst")

    # 6. Figures
    print("\n[INFO] Generating figures...")
    fig1_real_spreads(spread_panel, event_df)
    fig2_avalanche_powerlaw(event_df)
    fig3_networks(nets_corr, nets_partial, nets_mst)
    fig4_geometry([(snap_corr,"corr"),(snap_partial,"partial"),(snap_mst,"mst")])
    fig5_alpha_iet(soc_iet)
    fig6_alpha_kappa(soc_iet, country_corr)
    fig7_avalanche_timeline(event_df)
    fig8_dspread_diagnostic(spread_panel)

    print("\n" + "=" * 80)
    print("PIPELINE V3.0 COMPLETE.")
    print(f"Figures: {FIG_DIR}")
    print(f"Data:    {DATA_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    main()
