"""
Paper Figure 8: Monthly Distribution of LSTID-Active and Non-LSTID Days, 2015-2022

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
Bar chart of the number of LSTID-active (red) and non-LSTID (blue) days in each calendar month,
aggregated over 2015-2022. LSTID-active: amplitude >= 33 km with the QC below.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv

OUTPUT (written next to this script)
  20150101-20221231_monthly_event_counts_combined.png

RUN
  python fig8/fig8_monthly_lstid_counts.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 5; logic unchanged
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
# PUBLICATION FIGURE: Monthly Event Counts (All Years Combined) - SINGLE PANEL
# One panel: LSTID-active (red) and non-LSTID (blue) bars side by side per month
# (internal label "Calm" = non-LSTID day)
# N values shown as text above each bar
# ======================================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# SETTINGS
# =========================
CSV_FILE = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
OUTPUT_DIR = HERE

START_DATE = "2015-01-01"
END_DATE   = "2022-12-31"

R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5

AMPLITUDE_COL = "amplitude_km"
TID_THRESHOLD = 33.0

SHOW_FIGURE = True
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================================================================================
# LOAD + QC
# ======================================================================================
print("Loading base amplitude + applying QC...")

df = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()
df.index = pd.to_datetime(df.index, errors="raise")

date_min = pd.to_datetime(START_DATE)
date_max = pd.to_datetime(END_DATE)
df = df.loc[(df.index >= date_min) & (df.index <= date_max)].copy()

# ALWAYS filter low-R2 values
df[R2_COL] = pd.to_numeric(df[R2_COL], errors="raise")
df = df.loc[df[R2_COL] >= R2_MIN].copy()
df = df.loc[pd.to_numeric(df["T_hr"]).between(PERIOD_MIN, PERIOD_MAX, inclusive="both")].copy()

df[AMPLITUDE_COL] = pd.to_numeric(df[AMPLITUDE_COL], errors="raise")
df = df.loc[df[AMPLITUDE_COL].notna()].copy()

df["event_class"] = pd.Series(index=df.index, dtype="object")
df.loc[df[AMPLITUDE_COL] < TID_THRESHOLD, "event_class"] = "Calm"
df.loc[df[AMPLITUDE_COL] >= TID_THRESHOLD, "event_class"] = "TID"

print(
    f"After QC: total={len(df)}, "
    f"TID={(df['event_class'] == 'TID').sum()}, "
    f"Calm={(df['event_class'] == 'Calm').sum()}"
)

# ======================================================================================
# MONTHLY COUNTS
# ======================================================================================
df["month"] = df.index.month

tid_counts = df[df["event_class"] == "TID"].groupby("month").size()
calm_counts = df[df["event_class"] == "Calm"].groupby("month").size()

months = pd.Index(range(1, 13), name="month")
tid_counts = tid_counts.reindex(months, fill_value=0)
calm_counts = calm_counts.reindex(months, fill_value=0)

month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
x = np.arange(1, 13)

# ======================================================================================
# PLOT - single panel, grouped bars
# ======================================================================================
print("Creating figure...")

bar_width = 0.38
offsets = [-bar_width / 2, bar_width / 2]

fig, ax = plt.subplots(figsize=(16, 8))

bars_tid = ax.bar(
    x + offsets[0],
    tid_counts.values,
    width=bar_width,
    color="red",
    alpha=0.80,
    edgecolor="black",
    label=f"LSTID (N={int(tid_counts.sum())})"
)

bars_calm = ax.bar(
    x + offsets[1],
    calm_counts.values,
    width=bar_width,
    color="blue",
    alpha=0.80,
    edgecolor="black",
    label=f"Non-LSTID (N={int(calm_counts.sum())})"
)

# N labels above each bar
for bar in bars_tid:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + 0.5,
        str(int(h)),
        ha="center",
        va="bottom",
        fontsize=10,
        color="darkred",
        fontweight="bold"
    )

for bar in bars_calm:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + 0.5,
        str(int(h)),
        ha="center",
        va="bottom",
        fontsize=10,
        color="darkblue",
        fontweight="bold"
    )

ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=16)
ax.set_ylabel("Count (Days)", fontsize=18, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.tick_params(axis="y", labelsize=16)
ax.legend(fontsize=15, loc="upper right")

ax.set_title(
    "Monthly Distribution of LSTID-Active and Non-LSTID Days\n"
    f"{date_min.strftime('%Y-%m-%d')} to {date_max.strftime('%Y-%m-%d')}",
    fontsize=22,
    fontweight="bold",
    loc="center",
    pad=12
)

plt.tight_layout()

# Save
sDate_str = date_min.strftime("%Y%m%d")
eDate_str = date_max.strftime("%Y%m%d")
png_fpath = os.path.join(
    OUTPUT_DIR,
    f"{sDate_str}-{eDate_str}_monthly_event_counts_combined.png"
)

print(f"Saving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")

if SHOW_FIGURE:
    pass  # (notebook displayed the figure here)
print("Done!")
