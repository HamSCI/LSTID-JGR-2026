"""
Paper Figure 10: Correlation Analysis Between Minimum Skip Distance and Solar Flux (F10.7)

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
(a) daily minimum skip distance and F10.7, (b) daily scatter with Pearson r,
(c) both series after a 365-day centered running mean (removes the annual cycle),
(d) scatter of the 365-day means colored by year, with r and a 95% block-bootstrap CI
    (block length 365 d = smoothing window, 2000 resamples, seed 42; only ~7 independent samples),
(e) seasonal residual (daily value minus its 365-day mean) by day of year, its day-of-year mean
    (+-15 d window), and a least-squares linear fit to the cosine of the noon solar zenith angle at
    39.8N 98.6W, with the fit R^2 and the steepest 10-day changes.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  omni2_daily_3hWFoJxMID.lst
  sza_daily.csv

OUTPUT (written next to this script)
  fig_correlation_panel_minskip_flux.png

RUN
  python fig10/fig10_minskip_vs_f107.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 1; logic unchanged
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
# 5-PANEL PUBLICATION FIGURE: MIN SKIP DISTANCE VS SOLAR FLUX (F10.7)
#
# PANELS:
#   a) Daily time series
#   b) Daily correlation (native cadence)
#   c) Daily + low-pass mean time series
#   d) Low-pass means colored by year
#   e) Seasonal residual (daily minus 365-day mean) vs noon cos(SZA) fit
#
# OUTPUT:
#   fig_correlation_panel_minskip_flux.png
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
OMNI_LST   = os.path.join(DATA_DIR, "omni2_daily_3hWFoJxMID.lst")
SZA_CSV    = os.path.join(DATA_DIR, "sza_daily.csv")      # noon SZA at 39.83N 98.58W (gen_sza.py)
CLIM_HALF_WINDOW = 15                # +-days for the day-of-year climatology
SLOPE_LAG_DAYS   = 10
OUTPUT_DIR = HERE

START_DATE = "2015-01-01"
END_DATE   = "2022-12-31"

R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5

F107_FILL = 999.9   # OMNI fill value for missing F10.7

LOWPASS_WINDOW_FLUX = 365   # V11: 365 d nulls the annual cycle (180 d passed ~65% of it)
BLOCK_LENGTH_FLUX   = 365   # matched to LOWPASS_WINDOW_FLUX, as the paper states (was 30)
N_BOOTSTRAP_FLUX    = 2000
BOOT_SEED_FLUX      = 42

TICK_LABELSIZE = 16

RAW_ALPHA_ROW1_BLUE = 0.65
RAW_ALPHA_ROW1_RED  = 0.85
RAW_LW_ROW1_BLUE    = 0.8
RAW_LW_ROW1_RED     = 2.0

RAW_ALPHA_BG = 0.10
RAW_LW_BG    = 0.5

SMOOTH_LW    = 3.0
SMOOTH_ALPHA = 0.90

FIGSIZE = (20, 32)

# =========================
# PANEL-SPECIFIC AXIS LIMITS
# =========================

# Panel a
MIN_SKIP_YMIN_A = 200.0
MIN_SKIP_YMAX_A = 1600.0
F107_YMIN_A     = 50.0
F107_YMAX_A     = 200.0

# Panel b
MIN_SKIP_YMIN_B = 400.0
MIN_SKIP_YMAX_B = 1600.0
F107_XMIN_B     = 60.0
F107_XMAX_B     = 200.0

# Panel c
MIN_SKIP_YMIN_C = 400.0
MIN_SKIP_YMAX_C = 1600.0
F107_YMIN_C     = 50.0
F107_YMAX_C     = 200.0

# Panel d
MIN_SKIP_YMIN_D = 600.0
MIN_SKIP_YMAX_D = 1400.0
F107_XMIN_D     = 65.0
F107_XMAX_D     = 140.0

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
# LOAD F10.7 (DAILY)
# ======================================================================================
df_f107 = pd.read_csv(
    OMNI_LST,
    sep=r"\s+",
    header=None,
    names=["year", "day", "placeholder", "proton_density", "f10.7_index"]
)

df_f107 = df_f107.drop(columns=["placeholder"])
df_f107["date"] = pd.to_datetime(
    df_f107["year"].astype(str) + df_f107["day"].astype(str),
    format="%Y%j"
)
df_f107 = df_f107.set_index("date").sort_index()

df_f107["f10.7_index"] = pd.to_numeric(df_f107["f10.7_index"])
df_f107.loc[
    df_f107["f10.7_index"] >= F107_FILL,
    "f10.7_index"
] = np.nan

# ======================================================================================
# MERGE + DAILY GRID
# ======================================================================================
df_flux_full = df_base.merge(df_f107[["f10.7_index"]], left_index=True, right_index=True, how="inner")
df_flux_full = df_flux_full.sort_index().asfreq("D")

print(f"Flux daily grid: {len(df_flux_full)} days from {df_flux_full.index.min()} to {df_flux_full.index.max()}")

# ======================================================================================
# LOW-PASS ROLLING MEANS
# ======================================================================================
min_periods_flux = max(1, LOWPASS_WINDOW_FLUX // 2)

df_flux_full["ms_smooth"] = df_flux_full["min_combined_fit"].rolling(
    LOWPASS_WINDOW_FLUX, center=True, min_periods=min_periods_flux
).mean()

df_flux_full["flux_smooth"] = df_flux_full["f10.7_index"].rolling(
    LOWPASS_WINDOW_FLUX, center=True, min_periods=min_periods_flux
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

flux_daily_pairs = df_flux_full[["min_combined_fit", "f10.7_index"]].dropna()
flux_raw_idx = flux_daily_pairs.index
r_flux_daily = pearsonr(flux_daily_pairs["min_combined_fit"], flux_daily_pairs["f10.7_index"])[0]
n_flux_daily = len(flux_daily_pairs)

flux_smooth_pairs = df_flux_full.loc[flux_raw_idx, ["ms_smooth", "flux_smooth"]].dropna()

r_flux_smooth, flux_ci_low, flux_ci_high, n_flux_smooth = block_bootstrap_correlation(
    flux_smooth_pairs["ms_smooth"],
    flux_smooth_pairs["flux_smooth"],
    block_length=BLOCK_LENGTH_FLUX,
    n_bootstrap=N_BOOTSTRAP_FLUX,
    seed=BOOT_SEED_FLUX
)

print(f"Flux daily:  r = {r_flux_daily:.3f}, N = {n_flux_daily}")
print(f"Flux smooth: r = {r_flux_smooth:.3f}, 95% CI = [{flux_ci_low:.3f}, {flux_ci_high:.3f}], N = {n_flux_smooth}")

# ======================================================================================
# SHARED PREP FOR PANEL d
# ======================================================================================
x_all = flux_smooth_pairs["flux_smooth"].to_numpy()
y_all = flux_smooth_pairs["ms_smooth"].to_numpy()

z_all = np.polyfit(x_all, y_all, 1)
x_line_all = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
y_line_all = z_all[0] * x_line_all + z_all[1]

flux_smooth_pairs_year = flux_smooth_pairs.copy()
flux_smooth_pairs_year["year"] = flux_smooth_pairs_year.index.year
years = np.sort(flux_smooth_pairs_year["year"].unique())

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
gs = fig.add_gridspec(5, 1, hspace=0.42)
fig.subplots_adjust(top=0.93)

# ======================================================================================
# a) Daily time series
# ======================================================================================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.plot(
    df_flux_full.index, df_flux_full["min_combined_fit"],
    color="blue",
    alpha=RAW_ALPHA_ROW1_BLUE,
    linewidth=RAW_LW_ROW1_BLUE,
    linestyle="--",
    label="Daily"
)
ax_a.set_ylabel("Min Skip Distance [km]", fontdict={"weight": "bold", "size": 20}, color="blue")
ax_a.tick_params(axis="y", labelcolor="blue", labelsize=TICK_LABELSIZE)
ax_a.set_ylim(MIN_SKIP_YMIN_A, MIN_SKIP_YMAX_A)
style_time_axis(ax_a)
ax_a.legend(loc="upper left", fontsize=14)
ax_a.set_title("a) Daily Time Series", fontsize=20, fontweight="bold", loc="left", pad=10)

ax_a_t = ax_a.twinx()
ax_a_t.plot(
    df_flux_full.index, df_flux_full["f10.7_index"],
    color="red",
    alpha=RAW_ALPHA_ROW1_RED,
    linewidth=RAW_LW_ROW1_RED,
    linestyle="--",
    zorder=5,
    label="Daily"
)
ax_a_t.set_ylabel("Solar Flux (F10.7)", fontdict={"weight": "bold", "size": 20}, color="red")
ax_a_t.set_ylim(F107_YMIN_A, F107_YMAX_A)
emphasize_right_axis(ax_a_t, color="red")
ax_a_t.legend(loc="upper right", fontsize=14)

# ======================================================================================
# b) Daily correlation
# ======================================================================================
ax_b = fig.add_subplot(gs[1, 0])
ax_b.scatter(
    flux_daily_pairs["f10.7_index"], flux_daily_pairs["min_combined_fit"],
    alpha=0.30, s=5, color="navy"
)
ax_b.set_xlabel("Solar Flux (F10.7)", fontsize=18, fontweight="bold")
ax_b.set_ylabel("Min Skip Distance [km]", fontsize=18, fontweight="bold")
ax_b.set_xlim(F107_XMIN_B, F107_XMAX_B)
ax_b.set_ylim(MIN_SKIP_YMIN_B, MIN_SKIP_YMAX_B)
ax_b.grid(True, alpha=0.3)
ax_b.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
ax_b.text(
    0.95, 0.95,
    f"r = {r_flux_daily:.3f}\nN = {n_flux_daily} days",
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
    df_flux_full.index, df_flux_full["min_combined_fit"],
    color="blue", alpha=RAW_ALPHA_BG, linewidth=RAW_LW_BG, label="Daily"
)
ax_c.plot(
    df_flux_full.index, df_flux_full["ms_smooth"],
    color="blue", lw=SMOOTH_LW, alpha=SMOOTH_ALPHA, label=f"{LOWPASS_WINDOW_FLUX}-day mean"
)
ax_c.set_ylabel("Min Skip Distance [km]", fontdict={"weight": "bold", "size": 20}, color="blue")
ax_c.tick_params(axis="y", labelcolor="blue", labelsize=TICK_LABELSIZE)
ax_c.set_ylim(MIN_SKIP_YMIN_C, MIN_SKIP_YMAX_C)
style_time_axis(ax_c)
ax_c.legend(loc="upper left", fontsize=14)
ax_c.set_title(f"c) {LOWPASS_WINDOW_FLUX}-Day Mean Time Series", fontsize=20, fontweight="bold", loc="left", pad=10)

ax_c_t = ax_c.twinx()
ax_c_t.plot(
    df_flux_full.index, df_flux_full["f10.7_index"],
    color="red", alpha=RAW_ALPHA_BG, linewidth=RAW_LW_BG, label="Daily"
)
ax_c_t.plot(
    df_flux_full.index, df_flux_full["flux_smooth"],
    color="red", lw=SMOOTH_LW, alpha=SMOOTH_ALPHA, linestyle="--", label=f"{LOWPASS_WINDOW_FLUX}-day mean"
)
ax_c_t.set_ylabel("Solar Flux (F10.7)", fontdict={"weight": "bold", "size": 20}, color="red")
ax_c_t.set_ylim(F107_YMIN_C, F107_YMAX_C)
ax_c_t.tick_params(axis="y", labelcolor="red", labelsize=TICK_LABELSIZE)
ax_c_t.legend(loc="upper right", fontsize=14)

# ======================================================================================
# d) Low-pass means colored by year
# ======================================================================================
ax_d = fig.add_subplot(gs[3, 0])

cmap = plt.cm.get_cmap("viridis", len(years))

for i, yr in enumerate(years):
    dfy = flux_smooth_pairs_year.loc[flux_smooth_pairs_year["year"] == yr]
    ax_d.scatter(
        dfy["flux_smooth"], dfy["ms_smooth"],
        s=18, alpha=0.65, color=cmap(i), label=str(yr),
        edgecolors="none"
    )


ax_d.plot(x_line_all, y_line_all, "r--", lw=2, alpha=0.90)

ax_d.set_xlabel("Solar Flux (F10.7)", fontsize=18, fontweight="bold")
ax_d.set_ylabel("Min Skip Distance [km]", fontsize=18, fontweight="bold")
ax_d.set_xlim(F107_XMIN_D, F107_XMAX_D)
ax_d.set_ylim(MIN_SKIP_YMIN_D, MIN_SKIP_YMAX_D)
ax_d.grid(True, alpha=0.3)
ax_d.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)

text_d = (
    f"r = {r_flux_smooth:.3f}\n"
    f"95% CI = [{flux_ci_low:.3f}, {flux_ci_high:.3f}]\n"
    f"N = {n_flux_smooth}"
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
    f"d) {LOWPASS_WINDOW_FLUX}-Day Low-Pass Means Colored by Year",
    fontsize=20, fontweight="bold", loc="left", pad=10
)

# ======================================================================================
# e) Seasonal residual (daily skip minus its 365 d mean) vs a noon cos(SZA) fit
# ======================================================================================
df_sza = pd.read_csv(SZA_CSV, parse_dates=["date"], index_col="date")
df_flux_full["cos_sza"] = np.cos(np.deg2rad(df_sza["sza_noon_deg"])).reindex(df_flux_full.index)
df_flux_full["ms_resid"] = df_flux_full["min_combined_fit"] - df_flux_full["ms_smooth"]

doy_all = df_flux_full.index.dayofyear.to_numpy().copy()
doy_all[df_flux_full.index.is_leap_year & (doy_all > 59)] -= 1
doy_all = np.clip(doy_all, 1, 365)

def doy_climatology(values):
    out = np.full(365, np.nan)
    for k in range(1, 366):
        dist = np.minimum(np.abs(doy_all - k), 365 - np.abs(doy_all - k))
        m = (dist <= CLIM_HALF_WINDOW) & ~np.isnan(values)
        out[k - 1] = values[m].mean()
    return out

clim_resid = doy_climatology(df_flux_full["ms_resid"].to_numpy())
clim_cos   = doy_climatology(df_flux_full["cos_sza"].to_numpy())
fit_coef   = np.polyfit(clim_cos, clim_resid, 1)
clim_fit   = np.polyval(fit_coef, clim_cos)
r2_fit     = 1 - np.var(clim_resid - clim_fit) / np.var(clim_resid)

def max_abs_slope(c):
    cc = np.concatenate([c, c[:SLOPE_LAG_DAYS]])
    return np.max(np.abs(cc[SLOPE_LAG_DAYS:] - cc[:-SLOPE_LAG_DAYS]))

slope_obs = max_abs_slope(clim_resid)
slope_fit = max_abs_slope(clim_fit)

print(f"Seasonal residual: cos(SZA) fit R2 = {r2_fit:.2f}; max |slope| per {SLOPE_LAG_DAYS} d: observed {slope_obs:.0f} km, fit {slope_fit:.0f} km")
month_starts = pd.date_range("2021-01-01", "2021-12-01", freq="MS").dayofyear.to_numpy()
mon_obs = [np.nanmean(clim_resid[month_starts[i] - 1:(month_starts[i + 1] - 1 if i < 11 else 365)]) for i in range(12)]
mon_fit = [np.nanmean(clim_fit[month_starts[i] - 1:(month_starts[i + 1] - 1 if i < 11 else 365)]) for i in range(12)]
print("Monthly residual climatology, observed minus fit [km]:", [round(o - f) for o, f in zip(mon_obs, mon_fit)])

ax_e = fig.add_subplot(gs[4, 0])
doy_axis = np.arange(1, 366)
ax_e.scatter(doy_all, df_flux_full["ms_resid"], s=4, color="blue", alpha=0.12, label="Daily residual")
ax_e.plot(doy_axis, clim_resid, color="blue", lw=SMOOTH_LW, label=f"Day-of-year mean (±{CLIM_HALF_WINDOW} d)")
ax_e.plot(doy_axis, clim_fit, color="black", lw=2.5, linestyle="--", label="Linear fit to noon cos(SZA)")
ax_e.axhline(0, color="gray", lw=1)
ax_e.set_xlim(1, 365)
ax_e.set_ylim(-600, 500)
ax_e.set_xticks(month_starts)
ax_e.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
ax_e.set_ylabel("Skip Residual [km]", fontsize=18, fontweight="bold")
ax_e.set_xlabel("Day of Year", fontsize=18, fontweight="bold")
ax_e.grid(True, alpha=0.3)
ax_e.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
ax_e.legend(loc="lower left", fontsize=14, markerscale=4.0)
ax_e.text(
    0.98, 0.05,
    f"cos(SZA) fit R² = {r2_fit:.2f}\n"
    f"max slope: observed {slope_obs:.0f} km / {SLOPE_LAG_DAYS} d\n"
    f"max slope: fit {slope_fit:.0f} km / {SLOPE_LAG_DAYS} d",
    transform=ax_e.transAxes, fontsize=15,
    verticalalignment="bottom", horizontalalignment="right",
    bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.85)
)
ax_e.set_title(
    f"e) Seasonal Residual (Daily minus {LOWPASS_WINDOW_FLUX}-Day Mean) by Day of Year",
    fontsize=20, fontweight="bold", loc="left", pad=10
)

# ======================================================================================
# MAIN TITLE + SAVE
# ======================================================================================
fig.suptitle(
    "Correlation Analysis Between Minimum Skip Distance and Solar Flux (F10.7)",
    fontsize=26, fontweight="bold", y=0.965
)

png_fpath = os.path.join(OUTPUT_DIR, "fig_correlation_panel_minskip_flux.png")
print(f"Saving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")

pass  # (notebook displayed the figure here)
plt.close(fig)

print("\nDone.")
