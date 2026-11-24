"""
Paper Figure 6: Monthly Climatology of LSTID Amplitude and Candidate Drivers, 2015-2022

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
Monthly means of (a) fitted amplitude, (b) polar vortex index, (c) AE daily mean,
(d) SYM-H daily mean. Thin lines are individual years, the bold line is the all-year mean, and the
red dashed line marks the most active month (maximum; minimum for SYM-H). The polar vortex index
ends in July 2022 and AE is unavailable (fill values) for 2020.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  merra2_pvIndex.csv
  sym_ae.lst

OUTPUT (written next to this script)
  fig_seasonal_phase_amp_pv_ae_symh.png

RUN
  python fig6/fig6_seasonal_phase.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 7; logic unchanged
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
# SEASONAL PHASE FIGURE (V11, new): monthly climatology 2015-2022
#
# PANELS:
#   a) Fitted LSTID amplitude
#   b) Polar vortex index
#   c) AE index (daily mean)
#   d) SYM-H index (daily mean)
# Thin lines = individual years, bold line = all-year monthly mean.
# Red dashed line = most active month (maximum, or minimum for SYM-H, which goes negative with activity).
# Supports the phase argument: vortex in phase with the winter amplitude maximum,
# AE (canonical LSTID driver) peaks near equinox and is anti-phased with it.
#
# OUTPUT:
#   fig_seasonal_phase_amp_pv_ae_symh.png
#
# HARD FAIL RULE:
# - No guards, no fallbacks, no warnings. If columns or files are wrong, it crashes.
# ======================================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# SETTINGS
# =========================
CSV_FILE   = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
PV_CSV     = os.path.join(DATA_DIR, "merra2_pvIndex.csv")
SYM_AE_LST = os.path.join(DATA_DIR, "sym_ae.lst")
OUTPUT_DIR = HERE

START_DATE = "2015-01-01"
END_DATE   = "2022-12-31"

R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5

PV_MIN = -3.0
PV_MAX = 3.0
SYM_AE_FILL = 99999

TICK_LABELSIZE = 16
FIGSIZE = (14, 20)
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================================================================================
# LOAD + QC
# ======================================================================================
df = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()
for col in ["amplitude_km", "T_hr", R2_COL]:
    df[col] = pd.to_numeric(df[col])
df = df.loc[START_DATE:END_DATE]
df = df.loc[df[R2_COL] >= R2_MIN]
df = df.loc[df["T_hr"].between(PERIOD_MIN, PERIOD_MAX, inclusive="both")]
amp = df["amplitude_km"]

pv = pd.read_csv(PV_CSV, header=None, names=["date", "pv"])
pv["date"] = pd.to_datetime(pv["date"], format="%m/%d/%Y")
pv = pv.set_index("date").sort_index()["pv"].loc[START_DATE:END_DATE]
pv = pv.where((pv >= PV_MIN) & (pv <= PV_MAX))

sa = pd.read_csv(SYM_AE_LST, sep=r"\s+", comment="#", header=None, names=["year", "doy", "hour", "minute", "ae", "symh"])
sa.index = (pd.to_datetime(sa["year"].astype(str) + sa["doy"].astype(str), format="%Y%j")
            + pd.to_timedelta(sa["hour"], unit="h") + pd.to_timedelta(sa["minute"], unit="m"))
sa = sa[["ae", "symh"]].replace(SYM_AE_FILL, np.nan)
sa_daily = sa.resample("D").mean().loc[START_DATE:END_DATE]

series = [
    ("a) Fitted LSTID Amplitude", "Amplitude [km]", amp, +1),
    ("b) Polar Vortex Index", "Polar Vortex Index", pv, +1),
    ("c) AE Index (Daily Mean)", "AE [nT]", sa_daily["ae"], +1),
    ("d) SYM-H Index (Daily Mean)", "SYM-H [nT]", sa_daily["symh"], -1),
]

# ======================================================================================
# FIGURE
# ======================================================================================
fig, axes = plt.subplots(4, 1, figsize=FIGSIZE, sharex=True)
fig.subplots_adjust(hspace=0.35, top=0.93)
years = np.arange(pd.Timestamp(START_DATE).year, pd.Timestamp(END_DATE).year + 1)
cmap = plt.get_cmap("viridis", len(years))
x = np.arange(1, 13)

print("Monthly climatology 2015-2022 (all-year mean):")
for ax, (title, ylabel, s, active_sign) in zip(axes, series):
    by_year = s.groupby([s.index.year, s.index.month]).mean().unstack()
    for i, yr in enumerate(years):
        if yr in by_year.index:
            ax.plot(x, by_year.loc[yr].reindex(x), color=cmap(i), lw=1.2, alpha=0.6, label=str(yr))
    clim = s.groupby(s.index.month).mean().reindex(x)
    ax.plot(x, clim, color="black", lw=3.5, marker="o", ms=8, label="2015–2022 mean")
    most_active = (active_sign * clim).idxmax()
    ax.axvline(most_active, color="red", lw=2, linestyle="--")
    ax.set_ylabel(ylabel, fontsize=18, fontweight="bold")
    ax.set_title(title, fontsize=20, fontweight="bold", loc="left", pad=8)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
    print(f"  {title:32s} max {MONTHS[clim.idxmax() - 1]}  min {MONTHS[clim.idxmin() - 1]}  most active {MONTHS[most_active - 1]}  "
          f"years per month min {int(by_year.notna().sum().min())}  values {np.round(clim.values, 2).tolist()}")

axes[-1].set_xticks(x)
axes[-1].set_xticklabels(MONTHS)
axes[-1].set_xlim(0.7, 12.3)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.91, 0.5), fontsize=14, frameon=True)
fig.suptitle("Monthly Climatology of LSTID Amplitude and Candidate Drivers, 2015–2022",
             fontsize=22, fontweight="bold", y=0.975)

png_fpath = os.path.join(OUTPUT_DIR, "fig_seasonal_phase_amp_pv_ae_symh.png")
print(f"Saving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")
pass  # (notebook displayed the figure here)
plt.close(fig)
