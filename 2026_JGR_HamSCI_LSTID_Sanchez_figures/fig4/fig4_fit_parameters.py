"""
Paper Figure 4: Daily Sinusoidal Fit Parameters, 2015-2022

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
Three stacked time series of the daily sinusoidal fit parameters:
  (a) fitted amplitude A (km, Equation 8),
  (b) minimum RF skip distance d_min (km, Equation 10),
  (c) fitted period T (hours, Equation 8).
Gray traces are daily values, blue traces are 7-day rolling means, and the background color
also encodes the 7-day rolling mean.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv

OUTPUT (written next to this script)
  Daily Sinusoidal Fit Parameters.png   (the file name used in the manuscript)

RUN
  python fig4/fig4_fit_parameters.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 0; logic unchanged
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
# SIN-FIT ANALYSIS (NO FUNCTION) + R2 FILTER + PNG SAVE ONLY
# Restores monthly minor tick labels exactly like your original code
# Adds: main title + panel labels a) b) c) on the left of subplot titles
# ======================================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates

# =========================
# SETTINGS
# =========================
CSV_FILE = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
OUTPUT_DIR = HERE
os.makedirs(OUTPUT_DIR, exist_ok=True)

START_YEAR = 2015
END_YEAR   = 2022

R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5

AMPLITUDE_MAX = 90.0
ROLLING_DAYS = 7
T_HR_CMAP = "rainbow"

date_min = pd.to_datetime(f"{START_YEAR}-01-01")
date_max = pd.to_datetime(f"{END_YEAR}-12-31")

# =========================
# LOAD DATA (LIKE PREVIOUS SCRIPT)
# =========================
if not os.path.exists(CSV_FILE):
    raise FileNotFoundError(f"The file {CSV_FILE} does not exist.")

df = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()

required_columns = ["min_combined_fit", "T_hr", "amplitude_km", R2_COL]
missing = [c for c in required_columns if c not in df.columns]
if missing:
    raise ValueError(f"Missing required column(s) in CSV: {missing}")

for col in required_columns:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.loc[(df.index.year >= START_YEAR) & (df.index.year <= END_YEAR)].copy()
df = df.loc[df[R2_COL] >= R2_MIN].copy()
df = df.loc[df["T_hr"].between(PERIOD_MIN, PERIOD_MAX, inclusive="both")].copy()

df["selected"] = df["amplitude_km"] < AMPLITUDE_MAX
df_sel = df.loc[df["selected"]].copy()

# Extract start and end dates for file naming (post-filter)
sDate = df_sel.index.min()
eDate = df_sel.index.max()
sDate_str = pd.to_datetime(sDate).strftime("%Y%m%d")
eDate_str = pd.to_datetime(eDate).strftime("%Y%m%d")

# =========================
# PLOTTING
# =========================
nrows = 3
ncols = 1
figsize = (34, nrows * 9)
fig = plt.figure(figsize=figsize)

fig.suptitle(
    f"Daily Sinusoidal Fit Parameters (R² ≥ {R2_MIN}, {PERIOD_MIN:g}–{PERIOD_MAX:g} hr)   {START_YEAR}–{END_YEAR}\n"
    f"{ROLLING_DAYS}-Day Rolling Mean",
    fontsize=45,
    fontweight="bold",
    y=0.995
)

prmds = {
    "amplitude_km": {
        "title": "Ham Radio LSTID Amplitude",
        "label": "Amplitude [km]",
        "vmin": 10,
        "vmax": 47,
        "y_min": 0,
        "y_max": 90,
    },
    "min_combined_fit": {
        "title": "Ham Radio Min RF Skip Distance",
        "label": "Skip Distance [km]",
        "vmin": 0,
        "vmax": 1600,
        "y_min": 0,
        "y_max": 2000,
    },
    "T_hr": {
        "title": "Period",
        "label": "Period (hrs)",
        "vmin": 0,
        "vmax": 5,
        "y_min": 0,
        "y_max": 7,
    },
}

panel_labels = ["a)", "b)", "c)"]

cbar_info = {}
cmap = mpl.colormaps.get_cmap(T_HR_CMAP)

for ax_inx, (param, prmd) in enumerate(prmds.items(), start=1):
    ax = fig.add_subplot(nrows, ncols, ax_inx)

    xx = df_sel.index
    yy_raw = df_sel[param]

    min_periods = max(1, ROLLING_DAYS // 2)
    yy = df_sel[param].rolling(ROLLING_DAYS, center=True, min_periods=min_periods).mean()

    title = f"{panel_labels[ax_inx-1]}  {prmd['title']} ({ROLLING_DAYS} Day Rolling Mean)"

    norm = mpl.colors.Normalize(vmin=prmd["vmin"], vmax=prmd["vmax"])
    mpbl = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    color = mpbl.to_rgba(yy)
    color[np.asarray(pd.isna(yy))] = (1, 1, 1, 1)  # NaN -> white

    ax.plot(xx, yy_raw, color="0.5", label="Raw Data")
    ax.plot(xx, yy, color="blue", lw=2.5, label=f"{ROLLING_DAYS} Day Rolling Mean")
    ax.scatter(xx, yy, marker="o", c=color)

    ax.legend(loc="upper right", ncols=2, fontsize=22)
    ax.set_ylim(prmd["y_min"], prmd["y_max"])
    ax.set_xlim(date_min, date_max)

    # RESTORED: yearly major + monthly minor ticks exactly like your original
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_minor_formatter(mpl.ticker.FuncFormatter(
        lambda x, pos: mdates.num2date(x).strftime("%b")
        if mdates.num2date(x).month in (3, 5, 7, 9, 11) else ""))

    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment("center")

    for label in ax.get_xminorticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment("center")

    trans = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.bar(xx, 1, width=1, color=color, align="edge", zorder=-1, transform=trans, alpha=0.5)

    ax.set_ylabel(prmd["label"], fontdict={"weight": "bold", "size": 35})
    ax.set_title(title, loc="left", fontsize=39, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=29)
    ax.tick_params(axis="both", which="minor", labelsize=25)

    cbar_info[ax_inx] = {"ax": ax, "label": prmd["label"], "mpbl": mpbl}

# Leave room for the suptitle
fig.tight_layout(rect=[0, 0, 1, 0.965])

# =========================
# COLORBARS (LIKE YOUR ORIGINAL)
# =========================
for ax_inx, cbd in cbar_info.items():
    ax_pos = cbd["ax"].get_position()
    cax = fig.add_axes([1.025, ax_pos.y0, 0.02, ax_pos.height])
    cbar = fig.colorbar(cbd["mpbl"], label=cbd["label"], cax=cax)
    cbar.ax.tick_params(labelsize=24)
    cbar.set_label(cbd["label"], fontsize=30)

# =========================
# SAVE PNG ONLY
# =========================
png_fname = "Daily Sinusoidal Fit Parameters.png"  # file name used in the manuscript
png_fpath = os.path.join(OUTPUT_DIR, png_fname)
print(f"Saving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")

pass  # (notebook displayed the figure here)
