"""
Paper Figure 11: Minimum Skip Distance versus Sporadic-E Occurrence and foF2, by Season

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
2x2 panels of the seasonal component of minimum skip distance against
(a) sporadic-E occurrence in summer (May-Aug), (b) sporadic-E occurrence in winter (Nov-Feb),
(c) foF2 in summer, (d) foF2 in winter.
Sporadic-E occurrence: percent of daytime (13-23 UTC) ionograms from six GIRO stations with
foEs >= 4 MHz (missing foEs counted as no occurrence). foF2: daily mean over the same stations
and window. All series: 540-day running mean removed, then 30-day running mean; 95% block-bootstrap
CI with 30-day blocks (= smoothing window).

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  giro/giro_<STATION>_<YEAR>.txt (BC840, AU930, EG931, AL945, WP937, MHJ45)

OUTPUT (written next to this script)
  fig_correlation_panel_minskip_fof2_es_summer_winter.png

RUN
  python fig11/fig11_minskip_vs_ionosonde.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 2; logic unchanged
except file paths and headless plotting. Printed statistics match the values in the paper.
"""
import os
import matplotlib
matplotlib.use("Agg")   # headless: write files, no display

HERE = os.path.dirname(os.path.abspath(__file__))
# Input data folder: set LSTID_DATA_DIR to where the data lives (see data/README.md);
# defaults to ../data/ next to this script.
DATA_DIR = os.environ.get("LSTID_DATA_DIR", os.path.join(HERE, "..", "data"))

# ======================================================================================
# 2x2 SUMMARY FIGURE: WHICH REFLECTOR CONTROLS MINIMUM SKIP DISTANCE, BY SEASON
#
# Minimum skip distance is set by the densest reflecting layer available to return the
# signal. The two candidates at 14 MHz are the F2 layer and sporadic-E. This figure
# tests both against both seasons:
#
#            |  Summer (May-Aug)  |  Winter (Nov-Feb)
#   ---------+--------------------+-------------------
#   Es       |   (a)              |   (b)
#   foF2     |   (c)              |   (d)
#
# Filtering is IDENTICAL to the sporadic-E / F10.7 / polar-vortex correlation panels
# (1.5-yr high-pass + 30-day running mean, block-bootstrap 95% CI), computed on the
# FULL daily series; the season restriction is applied afterwards.
#
# OUTPUT: fig_correlation_panel_minskip_fof2_es_summer_winter.png (next to this script)
#
# HARD FAIL RULE: no guards, no fallbacks. Wrong columns/files -> crash.
# ======================================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# ============ SETTINGS (identical to the other correlation panels) ============
CSV_FILE   = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
OUTPUT_DIR = HERE
START_DATE, END_DATE = "2015-01-01", "2022-12-31"

R2_MIN = 0.35
CS_MIN, H0, H1 = 0, 13, 23
STATIONS  = ['BC840', 'AU930', 'EG931', 'AL945', 'WP937', 'MHJ45']
ES_THRESH = 4.0

HIGHPASS_DAYS, SMOOTH_DAYS = 540, 30
BLOCK_LENGTH, N_BOOTSTRAP, BOOT_SEED = 30, 2000, 42

SUMMER, WINTER = [5, 6, 7, 8], [11, 12, 1, 2]
TICK_LABELSIZE = 16

os.makedirs(OUTPUT_DIR, exist_ok=True)
date_min, date_max = pd.to_datetime(START_DATE), pd.to_datetime(END_DATE)

# ============ HamSCI daily min-skip ============
df = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()
df["min_combined_fit"] = pd.to_numeric(df["min_combined_fit"])
df["r2"] = pd.to_numeric(df["r2"])
df = df.loc[(df.index >= date_min) & (df.index <= date_max)]
df = df.loc[df["r2"] >= R2_MIN]
df = df.loc[pd.to_numeric(df["T_hr"]).between(1.0, 4.5, inclusive="both")]   # paper Sec 3.5 period QC
ms_daily = df["min_combined_fit"].resample("D").mean()

# ============ GIRO daily foEs occurrence + foF2 ============
rows = []
for stn in STATIONS:
    for yr in range(2015, 2023):
        with open(f"{DATA_DIR}/giro/giro_{stn}_{yr}.txt") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 24 or float(p[1]) < CS_MIN:
                    continue
                t = pd.Timestamp(p[0]).tz_convert(None)
                if t.hour < H0 or t.hour >= H1 or t < date_min or t > date_max:
                    continue
                foEs = np.nan if p[12] == "---" else float(p[12])
                foF2 = np.nan if p[2]  == "---" else float(p[2])
                rows.append((t.normalize(), foEs, foF2))
g = pd.DataFrame(rows, columns=["date", "foEs", "foF2"])
es_ser = g.groupby("date")["foEs"].apply(lambda s: (s.fillna(-9) >= ES_THRESH).sum() / len(s)) * 100.0
f2_ser = g.groupby("date")["foF2"].mean()

# ============ high-pass + smooth on the FULL series ============
grid = pd.date_range(START_DATE, END_DATE, freq="D")
highpass = lambda s, w: s - s.rolling(w, center=True, min_periods=max(1, w // 2)).mean()
smooth   = lambda s, w: s.rolling(w, center=True, min_periods=max(1, w // 2)).mean()

ms = ms_daily.reindex(grid)
ms_sm = smooth(highpass(ms, HIGHPASS_DAYS), SMOOTH_DAYS)
V_SM = {
    "es":   smooth(highpass(es_ser.reindex(grid), HIGHPASS_DAYS), SMOOTH_DAYS),
    "foF2": smooth(highpass(f2_ser.reindex(grid), HIGHPASS_DAYS), SMOOTH_DAYS),
}
RAW = {"es": es_ser.reindex(grid), "foF2": f2_ser.reindex(grid)}

def block_bootstrap_correlation(x, y, *, block_length=BLOCK_LENGTH,
                                n_bootstrap=N_BOOTSTRAP, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    c = pd.concat([x, y], axis=1).dropna()
    xa, ya = c.iloc[:, 0].to_numpy(), c.iloc[:, 1].to_numpy()
    n = len(c)
    r_obs = pearsonr(xa, ya)[0]
    starts = np.arange(0, max(1, n - block_length + 1))
    boot = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        idx = []
        while len(idx) < n:
            s = int(rng.choice(starts)); idx.extend(range(s, s + block_length))
        boot[b] = pearsonr(xa[np.array(idx[:n])], ya[np.array(idx[:n])])[0]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return r_obs, lo, hi, n

# ============ FIGURE ============
PANELS = [("es",   SUMMER, "a"), ("es",   WINTER, "b"),
          ("foF2", SUMMER, "c"), ("foF2", WINTER, "d")]
XLAB = {"es": "Sporadic-E Occurrence Residual [%]", "foF2": "foF2 Residual [MHz]"}
VNAME = {"es": "Sporadic-E", "foF2": "foF2"}
SNAME = {tuple(SUMMER): "Summer (May-Aug)", tuple(WINTER): "Winter (Nov-Feb)"}

ms_lim = np.nanpercentile(np.abs(ms_sm.dropna()), 99) * 1.15
xlims  = {v: np.nanpercentile(np.abs(V_SM[v].dropna()), 99) * 1.15 for v in ["es", "foF2"]}

years = np.arange(2015, 2023)
cmap = plt.cm.get_cmap("viridis", len(years))

fig, axes = plt.subplots(2, 2, figsize=(19, 15))
results = {}
for ax, (var, months, tag) in zip(axes.ravel(), PANELS):
    keep = pd.Series(grid.month.isin(months), index=grid)
    x, y = V_SM[var].where(keep), ms_sm.where(keep)
    r_obs, lo, hi, n = block_bootstrap_correlation(x, y)
    r_raw = pearsonr(*pd.concat([RAW[var].where(keep), ms.where(keep)],
                                axis=1).dropna().to_numpy().T)[0]
    results[(var, tag)] = (r_obs, lo, hi, n, r_raw)
    sig = not (lo < 0 < hi)

    d = pd.concat([x, y], axis=1, keys=["x", "y"]).dropna()
    d["year"] = d.index.year
    for i, yr in enumerate(years):
        dd = d[d["year"] == yr]
        ax.scatter(dd["x"], dd["y"], s=14, alpha=0.65, color=cmap(i),
                   label=str(yr) if tag == "a" else None, edgecolors="none")
    z = np.polyfit(d["x"], d["y"], 1)
    # clip the fit line to the span of the data so null panels do not imply
    # a relationship across regions where no observations exist
    xl = np.linspace(d["x"].min(), d["x"].max(), 200)
    ax.plot(xl, z[0] * xl + z[1], color="red" if sig else "grey",
            ls="--", lw=2.5 if sig else 1.8, alpha=0.9)

    ax.set_xlim(-xlims[var], xlims[var]); ax.set_ylim(-ms_lim, ms_lim)
    ax.axhline(0, color="k", lw=0.6, ls=":"); ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.grid(True, alpha=0.3); ax.tick_params(axis="both", labelsize=TICK_LABELSIZE)
    ax.set_xlabel(XLAB[var], fontsize=18, fontweight="bold")
    ax.set_ylabel("Min Skip Residual [km]", fontsize=18, fontweight="bold")
    ax.text(0.96, 0.96,
            f"r = {r_obs:.3f}\n95% CI [{lo:.2f}, {hi:.2f}]\n"
            f"N = {n}\nraw daily r = {r_raw:.3f}",
            transform=ax.transAxes, fontsize=15, va="top", ha="right",
            bbox=dict(boxstyle="round",
                      facecolor="lightblue" if sig else "lightgrey", alpha=0.9))
    ax.set_title(f"{tag}) {VNAME[var]}, {SNAME[tuple(months)]}",
                 fontsize=20, fontweight="bold", loc="left", pad=10)
    for sp in ax.spines.values():
        sp.set_linewidth(2.4 if sig else 0.8)
        sp.set_color("firebrick" if sig else "grey")

fig.suptitle("Minimum Skip Distance versus Sporadic-E Occurrence and foF2, by Season\n"
             f"(seasonal component (residual), {HIGHPASS_DAYS/365:.1f}-yr high-pass + {SMOOTH_DAYS}-day mean)",
             fontsize=24, fontweight="bold", y=0.978)
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, title="Year", fontsize=14, title_fontsize=14,
           loc="center right", bbox_to_anchor=(1.06, 0.5), frameon=True)
fig.tight_layout(rect=[0, 0, 1, 0.94])

png = os.path.join(OUTPUT_DIR, "fig_correlation_panel_minskip_fof2_es_summer_winter.png")
fig.savefig(png, dpi=300, bbox_inches="tight")
print("Saved:", png)
for (var, tag), (r_obs, lo, hi, n, r_raw) in results.items():
    print(f"  {tag}) {VNAME[var]:<11} r = {r_obs:+.3f}  CI [{lo:+.3f},{hi:+.3f}]  "
          f"N = {n:4d}   (raw daily r = {r_raw:+.3f})")
pass  # (notebook displayed the figure here)
plt.close(fig)
