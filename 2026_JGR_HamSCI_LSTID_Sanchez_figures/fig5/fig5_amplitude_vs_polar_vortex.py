"""
Paper Figure 5: Correlation Analysis Between Fitted Amplitude and the Polar Vortex Index

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
(a) daily time series of fitted amplitude and the polar vortex index,
(b) daily scatter with Pearson r,
(c) both series after a 60-day centered running mean,
(d) scatter of the 60-day means colored by year, with r and a 95% block-bootstrap CI
    (block length 60 d = smoothing window, 2000 resamples, seed 42).
A 60-day window retains the annual cycle, so (d) mainly measures agreement of two seasonal cycles.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  merra2_pvIndex.csv

OUTPUT (written next to this script)
  fig_correlation_panel_amp_pv.png

RUN
  python fig5/fig5_amplitude_vs_polar_vortex.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 3; logic unchanged
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
# 4-PANEL PUBLICATION FIGURE: FITTED AMPLITUDE VS POLAR VORTEX INDEX
#
# PANELS:
#   a) Daily time series
#   b) Daily correlation (native cadence)
#   c) Daily + low-pass mean time series
#   d) Low-pass means colored by year
#
# OUTPUT:
#   fig_correlation_panel_amp_pv.png
#
# HARD FAIL RULE:
# - No guards, no fallbacks, no warnings. If columns or files are wrong, it crashes.
# ======================================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import pearsonr

# =========================
# SETTINGS
# =========================
CSV_FILE   = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
PV_FILE    = os.path.join(DATA_DIR, "merra2_pvIndex.csv")
OUTPUT_DIR = HERE

START_DATE = "2015-01-01"
END_DATE   = "2022-12-31"

R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5

PV_MIN = -3.0
PV_MAX = 3.0

LOWPASS_WINDOW_PV = 60
BLOCK_LENGTH_PV   = 60
N_BOOTSTRAP_PV    = 2000
BOOT_SEED_PV      = 42

TICK_LABELSIZE = 16

RAW_ALPHA_ROW1_BLUE = 0.65
RAW_ALPHA_ROW1_RED  = 0.85
RAW_LW_ROW1_BLUE    = 0.8
RAW_LW_ROW1_RED     = 2.0

RAW_ALPHA_BG = 0.10
RAW_LW_BG    = 0.5

SMOOTH_LW    = 3.0
SMOOTH_ALPHA = 0.90

FIGSIZE = (20, 26)

# =========================
# PANEL-SPECIFIC AXIS LIMITS
# =========================

# Panel a
AMP_YMIN_A = -10.0
AMP_YMAX_A = 70.0
PV_YMIN_A  = -3.0
PV_YMAX_A  = 3.0

# Panel b
AMP_YMIN_B = 0.0
AMP_YMAX_B = 60.0
PV_XMIN_B  = -1.5
PV_XMAX_B  = 2.0

# Panel c
AMP_YMIN_C = 0.0
AMP_YMAX_C = 60.0
PV_YMIN_C  = -3.0
PV_YMAX_C  = 3.0

# Panel d
AMP_YMIN_D = 10.0
AMP_YMAX_D = 50.0
PV_XMIN_D  = -1.5
PV_XMAX_D  = 2.0

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================================================================================
# LOAD BASE + QC (STRICT)
# ======================================================================================
print("Loading data...")

df_base = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()

df_base["amplitude_km"] = pd.to_numeric(df_base["amplitude_km"])
df_base["min_combined_fit"] = pd.to_numeric(df_base["min_combined_fit"])
df_base[R2_COL] = pd.to_numeric(df_base[R2_COL])

date_min = pd.to_datetime(START_DATE)
date_max = pd.to_datetime(END_DATE)

df_base = df_base.loc[(df_base.index >= date_min) & (df_base.index <= date_max)].copy()
df_base = df_base.loc[df_base[R2_COL] >= R2_MIN].copy()
df_base["T_hr"] = pd.to_numeric(df_base["T_hr"])
df_base = df_base.loc[df_base["T_hr"].between(PERIOD_MIN, PERIOD_MAX, inclusive="both")].copy()

# ======================================================================================
# LOAD PV (DAILY)
# ======================================================================================
df_pv = pd.read_csv(PV_FILE, header=None, names=["date", "pv_index"])
df_pv["date"] = pd.to_datetime(df_pv["date"], format="%m/%d/%Y")
df_pv = df_pv.set_index("date").sort_index()

df_pv["pv_index"] = pd.to_numeric(df_pv["pv_index"])
df_pv.loc[
    (df_pv["pv_index"] < PV_MIN) |
    (df_pv["pv_index"] > PV_MAX),
    "pv_index"
] = np.nan

# ======================================================================================
# MERGE + DAILY GRID
# ======================================================================================
df_pv_full = df_base.merge(df_pv[["pv_index"]], left_index=True, right_index=True, how="inner")
df_pv_full = df_pv_full.sort_index().asfreq("D")

print(f"PV daily grid: {len(df_pv_full)} days from {df_pv_full.index.min()} to {df_pv_full.index.max()}")

# ======================================================================================
# LOW-PASS ROLLING MEANS
# ======================================================================================
min_periods_pv = max(1, LOWPASS_WINDOW_PV // 2)

df_pv_full["amp_smooth"] = df_pv_full["amplitude_km"].rolling(
    LOWPASS_WINDOW_PV, center=True, min_periods=min_periods_pv
).mean()

df_pv_full["pv_smooth"] = df_pv_full["pv_index"].rolling(
    LOWPASS_WINDOW_PV, center=True, min_periods=min_periods_pv
).mean()

# ======================================================================================
# BLOCK BOOTSTRAP (MOVING BLOCKS)
# ======================================================================================
def block_bootstrap_correlation(x: pd.Series, y: pd.Series, *, block_length: int, n_bootstrap: int, seed: int):
    rng = np.random.default_rng(seed)

    combined = pd.concat([x, y], axis=1).dropna()
    n = len(combined)

    x_arr = combined.iloc[:, 0].to_numpy()
    y_arr = combined.iloc[:, 1].to_numpy()

    r_obs = pearsonr(x_arr, y_arr)[0]

    starts = np.arange(0, n - block_length + 1)
    boot_rs = np.empty(n_bootstrap, dtype=float)

    for b in range(n_bootstrap):
        idx = []
        while len(idx) < n:
            s = int(rng.choice(starts))
            idx.extend(range(s, s + block_length))
        idx = np.array(idx[:n], dtype=int)
        boot_rs[b] = pearsonr(x_arr[idx], y_arr[idx])[0]

    ci_low, ci_high = np.nanpercentile(boot_rs, [2.5, 97.5])
    return r_obs, ci_low, ci_high, int(n)

# ======================================================================================
# CORRELATIONS
# ======================================================================================
print("\nComputing correlations...")

pv_daily_pairs = df_pv_full[["amplitude_km", "pv_index"]].dropna()
pv_raw_idx = pv_daily_pairs.index
r_pv_daily = pearsonr(pv_daily_pairs["amplitude_km"], pv_daily_pairs["pv_index"])[0]
n_pv_daily = len(pv_daily_pairs)

pv_smooth_pairs = df_pv_full.loc[pv_raw_idx, ["amp_smooth", "pv_smooth"]].dropna()

r_pv_smooth, pv_ci_low, pv_ci_high, n_pv_smooth = block_bootstrap_correlation(
    pv_smooth_pairs["amp_smooth"],
    pv_smooth_pairs["pv_smooth"],
    block_length=BLOCK_LENGTH_PV,
    n_bootstrap=N_BOOTSTRAP_PV,
    seed=BOOT_SEED_PV
)

print(f"PV daily:  r = {r_pv_daily:.3f}, N = {n_pv_daily}")
print(f"PV smooth: r = {r_pv_smooth:.3f}, 95% CI = [{pv_ci_low:.3f}, {pv_ci_high:.3f}], N = {n_pv_smooth}")

# ======================================================================================
# SHARED PREP FOR PANEL d
# ======================================================================================
x_all = pv_smooth_pairs["pv_smooth"].to_numpy()
y_all = pv_smooth_pairs["amp_smooth"].to_numpy()

z_all = np.polyfit(x_all, y_all, 1)
x_line_all = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
y_line_all = z_all[0] * x_line_all + z_all[1]

pv_smooth_pairs_year = pv_smooth_pairs.copy()
pv_smooth_pairs_year["year"] = pv_smooth_pairs_year.index.year
years = np.sort(pv_smooth_pairs_year["year"].unique())

# ======================================================================================
# FIGURE HELPERS
# ======================================================================================
def style_time_axis(ax):
    ax.set_xlim(date_min, date_max)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    for lab in ax.get_xticklabels(which="major"):
        lab.set_rotation(0)
        lab.set_horizontalalignment("center")

    ax.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
    ax.tick_params(axis="x", which="minor", labelsize=0)
    ax.grid(True, alpha=0.3)

def emphasize_right_axis(ax_t, color="red", spine_lw=2.6, tick_w=2.0):
    ax_t.spines["right"].set_color(color)
    ax_t.spines["right"].set_linewidth(spine_lw)
    ax_t.tick_params(axis="y", colors=color, width=tick_w, length=7, labelsize=TICK_LABELSIZE)

# ======================================================================================
# FIGURE
# ======================================================================================
print("\nCreating figure...")
fig = plt.figure(figsize=FIGSIZE)
gs = fig.add_gridspec(4, 1, hspace=0.42)
fig.subplots_adjust(top=0.91)

# ======================================================================================
# a) Daily time series
# ======================================================================================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.plot(
    df_pv_full.index, df_pv_full["amplitude_km"],
    color="blue",
    alpha=RAW_ALPHA_ROW1_BLUE,
    linewidth=RAW_LW_ROW1_BLUE,
    linestyle="--",
    label="Daily"
)
ax_a.set_ylabel("Fitted Amplitude [km]", fontdict={"weight": "bold", "size": 20}, color="blue")
ax_a.tick_params(axis="y", labelcolor="blue", labelsize=TICK_LABELSIZE)
ax_a.set_ylim(AMP_YMIN_A, AMP_YMAX_A)
style_time_axis(ax_a)
ax_a.legend(loc="upper left", fontsize=14)
ax_a.set_title("a) Daily Time Series", fontsize=20, fontweight="bold", loc="left", pad=10)

ax_a_t = ax_a.twinx()
ax_a_t.plot(
    df_pv_full.index, df_pv_full["pv_index"],
    color="red",
    alpha=RAW_ALPHA_ROW1_RED,
    linewidth=RAW_LW_ROW1_RED,
    linestyle="--",
    zorder=5,
    label="Daily"
)
ax_a_t.set_ylabel("Polar Vortex Index", fontdict={"weight": "bold", "size": 20}, color="red")
ax_a_t.set_ylim(PV_YMIN_A, PV_YMAX_A)
emphasize_right_axis(ax_a_t, color="red")
ax_a_t.legend(loc="upper right", fontsize=14)

# ======================================================================================
# b) Daily correlation
# ======================================================================================
ax_b = fig.add_subplot(gs[1, 0])
ax_b.scatter(
    pv_daily_pairs["pv_index"], pv_daily_pairs["amplitude_km"],
    alpha=0.30, s=5, color="navy"
)
ax_b.set_xlabel("Polar Vortex Index", fontsize=18, fontweight="bold")
ax_b.set_ylabel("Fitted Amplitude [km]", fontsize=18, fontweight="bold")
ax_b.set_xlim(PV_XMIN_B, PV_XMAX_B)
ax_b.set_ylim(AMP_YMIN_B, AMP_YMAX_B)
ax_b.grid(True, alpha=0.3)
ax_b.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
ax_b.text(
    0.95, 0.95,
    f"r = {r_pv_daily:.3f}\nN = {n_pv_daily} days",
    transform=ax_b.transAxes, fontsize=16,
    verticalalignment="top", horizontalalignment="right",
    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8)
)
ax_b.set_title("b) Daily Correlation (Native Cadence)", fontsize=20, fontweight="bold", loc="left", pad=10)

# ======================================================================================
# c) Daily + smoothed time series
# ======================================================================================
ax_c = fig.add_subplot(gs[2, 0])
ax_c.plot(
    df_pv_full.index, df_pv_full["amplitude_km"],
    color="blue", alpha=RAW_ALPHA_BG, linewidth=RAW_LW_BG, label="Daily"
)
ax_c.plot(
    df_pv_full.index, df_pv_full["amp_smooth"],
    color="blue", lw=SMOOTH_LW, alpha=SMOOTH_ALPHA, label=f"{LOWPASS_WINDOW_PV}-day mean"
)
ax_c.set_ylabel("Fitted Amplitude [km]", fontdict={"weight": "bold", "size": 20}, color="blue")
ax_c.tick_params(axis="y", labelcolor="blue", labelsize=TICK_LABELSIZE)
ax_c.set_ylim(AMP_YMIN_C, AMP_YMAX_C)
style_time_axis(ax_c)
ax_c.legend(loc="upper left", fontsize=14)
ax_c.set_title(f"c) {LOWPASS_WINDOW_PV}-Day Mean Time Series", fontsize=20, fontweight="bold", loc="left", pad=10)

ax_c_t = ax_c.twinx()
ax_c_t.plot(
    df_pv_full.index, df_pv_full["pv_index"],
    color="red", alpha=RAW_ALPHA_BG, linewidth=RAW_LW_BG, label="Daily"
)
ax_c_t.plot(
    df_pv_full.index, df_pv_full["pv_smooth"],
    color="red", lw=SMOOTH_LW, alpha=SMOOTH_ALPHA, linestyle="--", label=f"{LOWPASS_WINDOW_PV}-day mean"
)
ax_c_t.set_ylabel("Polar Vortex Index", fontdict={"weight": "bold", "size": 20}, color="red")
ax_c_t.set_ylim(PV_YMIN_C, PV_YMAX_C)
ax_c_t.tick_params(axis="y", labelcolor="red", labelsize=TICK_LABELSIZE)
ax_c_t.legend(loc="upper right", fontsize=14)

# ======================================================================================
# d) Low-pass means colored by year
# ======================================================================================
ax_d = fig.add_subplot(gs[3, 0])

cmap = plt.cm.get_cmap("viridis", len(years))

for i, yr in enumerate(years):
    dfy = pv_smooth_pairs_year.loc[pv_smooth_pairs_year["year"] == yr]
    ax_d.scatter(
        dfy["pv_smooth"], dfy["amp_smooth"],
        s=18, alpha=0.65, color=cmap(i), label=str(yr),
        edgecolors="none"
    )


ax_d.plot(x_line_all, y_line_all, "r--", lw=2, alpha=0.90)

ax_d.set_xlabel("Polar Vortex Index", fontsize=18, fontweight="bold")
ax_d.set_ylabel("Fitted Amplitude [km]", fontsize=18, fontweight="bold")
ax_d.set_xlim(PV_XMIN_D, PV_XMAX_D)
ax_d.set_ylim(AMP_YMIN_D, AMP_YMAX_D)
ax_d.grid(True, alpha=0.3)
ax_d.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)

text_d = (
    f"r = {r_pv_smooth:.3f}\n"
    f"95% CI = [{pv_ci_low:.3f}, {pv_ci_high:.3f}]\n"
    f"N = {n_pv_smooth}"
)

ax_d.text(
    0.95, 0.95,
    text_d,
    transform=ax_d.transAxes, fontsize=15,
    verticalalignment="top",
    horizontalalignment="right",
    bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.85)
)

ax_d.legend(
    title="Year",
    fontsize=14,
    title_fontsize=14,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=True,
    markerscale=3.0
)

ax_d.set_title(
    f"d) {LOWPASS_WINDOW_PV}-Day Low-Pass Means Colored by Year",
    fontsize=20, fontweight="bold", loc="left", pad=10
)

# ======================================================================================
# MAIN TITLE + SAVE
# ======================================================================================
fig.suptitle(
    "Correlation Analysis Between Fitted Amplitude and the Polar Vortex Index",
    fontsize=26, fontweight="bold", y=0.975
)

png_fpath = os.path.join(OUTPUT_DIR, "fig_correlation_panel_amp_pv.png")
print(f"Saving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")

pass  # (notebook displayed the figure here)
plt.close(fig)

print("\nDone.")
