"""
Paper Figure 7: Manual versus Automated Seasonal Validation for 2017

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
(a) RMSE between manual and automated monthly LSTID-day fractions for 2017 as a function of
    the automated amplitude threshold (10-90 km, 0.5 km steps); the minimum selects 33 km.
(b) 2017 monthly LSTID-day fractions, manual (red) vs automated at the selected threshold (gray).
Manual LSTID day: >= 1 hour of LSTID activity identified by inspection (column `tid`).
Automated LSTID day: amplitude >= threshold, period 1-4.5 hr, R^2 >= 0.35.
Also writes the full threshold sweep (with daily precision/recall/F1) as a CSV.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  2017_data_full_V2.csv

OUTPUT (written next to this script)
  fig_rmse_thresh.png

RUN
  python fig7/fig7_manual_validation_2017.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 4; logic unchanged
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
# VALIDATION FIGURE (2017): Monthly manual vs automated threshold sweep
#
# SINGLE FIGURE, 2 PANELS ONLY:
#   a) RMSE(raw monthly fractions) vs amplitude threshold
#   b) Raw monthly TID-day fractions for the selected threshold
#
# SIMPLE METHOD:
#   - Choose threshold using MONTHLY RMSE of RAW monthly TID-day fractions
#   - No mean-normalized optimization for threshold selection
#   - Still compute seasonal correlation and daily F1 as reference diagnostics
#
# Outputs:
#   1) One 2-panel PNG figure
#   2) CSV table of sweep metrics
#
# HARD FAIL RULE: no guards/fallbacks. If columns/files are wrong, it crashes.
# ======================================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import calendar

from sklearn.metrics import (
    f1_score,
    precision_score,
    recall_score,
    balanced_accuracy_score,
    accuracy_score,
    confusion_matrix,
)

# =========================
# SETTINGS
# =========================
# Manual 2017 (hours/day)
MANUAL_CSV = os.path.join(DATA_DIR, "2017_data_full_V2.csv")     # must contain: date, tid
MANUAL_DATE_COL = "date"
MANUAL_TID_HOURS_COL = "tid"
MANUAL_BINARY_HOUR_THRESHOLD = 1.0       # >= 1 hr => manual TID day

# Automated daily fit parameters
AUTO_CSV = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")  # daily index = date; must contain r2, amplitude_km, T_hr
AUTO_DATE_INDEX = True

YEAR = 2017

# Automated QC
R2_COL = "r2"
R2_MIN = 0.35

AMPLITUDE_COL = "amplitude_km"
MIN_AMPLITUDE = 0.0
MAX_AMPLITUDE = 90.0

PERIOD_COL = "T_hr"
PERIOD_MIN = 1.0
PERIOD_MAX = 4.5   # V11: match paper Sec 3.5 (was 5.0)

# -------------------------
# SWEEP SETTINGS
# -------------------------
THRESH_START = 10.0
THRESH_END   = 90.0
THRESH_STEP  = 0.5

OUTPUT_DIR = HERE
os.makedirs(OUTPUT_DIR, exist_ok=True)

SHOW_FIGURE = True

# -------------------------
# Plot styling
# -------------------------
FIGSIZE_2PANEL = (14, 16)
BAR_ALPHA = 0.80
EDGE_COLOR = "black"

TITLE_FONTSIZE = 20
SUPTITLE_FONTSIZE = 22
LABEL_FONTSIZE = 18
TICK_FONTSIZE = 16
BOX_FONTSIZE = 14
LEGEND_FONTSIZE = 14

month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
x = np.arange(1, 13)
months = pd.Index(range(1, 13), name="month")

# ======================================================================================
# LOAD MANUAL 2017
# ======================================================================================
dfm = pd.read_csv(MANUAL_CSV)
dfm[MANUAL_DATE_COL] = pd.to_datetime(dfm[MANUAL_DATE_COL], errors="raise")
dfm[MANUAL_TID_HOURS_COL] = pd.to_numeric(dfm[MANUAL_TID_HOURS_COL], errors="raise")

dfm = dfm.sort_values(MANUAL_DATE_COL).copy()
dfm = dfm.loc[dfm[MANUAL_DATE_COL].dt.year == YEAR].copy()
dfm["date"] = dfm[MANUAL_DATE_COL].dt.normalize()
dfm["month"] = dfm["date"].dt.month

# Manual binary day label
dfm["manual_tid_day"] = (dfm[MANUAL_TID_HOURS_COL] >= MANUAL_BINARY_HOUR_THRESHOLD).astype(int)

# Keep one row per day
manual_daily = (
    dfm[["date", "month", "manual_tid_day"]]
    .drop_duplicates(subset=["date"])
    .sort_values("date")
    .copy()
)

manual_tid_days = manual_daily.groupby("month")["manual_tid_day"].sum().reindex(months, fill_value=0)

manual_days_in_month = pd.Series(
    {m: calendar.monthrange(YEAR, m)[1] for m in range(1, 13)}
).reindex(months)

manual_frac = (manual_tid_days / manual_days_in_month).astype(float)
manual_mean = float(manual_frac.mean())
manual_norm = manual_frac / manual_mean

# ======================================================================================
# LOAD AUTOMATED 2017
# ======================================================================================
if AUTO_DATE_INDEX:
    dfa = pd.read_csv(AUTO_CSV, parse_dates=True, index_col=0).sort_index()
    dfa.index = pd.to_datetime(dfa.index, errors="raise")
    dfa["date"] = dfa.index.normalize()
else:
    dfa = pd.read_csv(AUTO_CSV)
    dfa["date"] = pd.to_datetime(dfa["date"], errors="raise").dt.normalize()
    dfa = dfa.sort_values("date").copy()

dfa = dfa.loc[dfa["date"].dt.year == YEAR].copy()

dfa[R2_COL] = pd.to_numeric(dfa[R2_COL], errors="raise")
dfa[AMPLITUDE_COL] = pd.to_numeric(dfa[AMPLITUDE_COL], errors="raise")
dfa[PERIOD_COL] = pd.to_numeric(dfa[PERIOD_COL], errors="raise")

# QC
dfa = dfa.loc[dfa[R2_COL] >= R2_MIN].copy()
dfa = dfa.loc[dfa[AMPLITUDE_COL].notna()].copy()
dfa = dfa.loc[(dfa[AMPLITUDE_COL] > MIN_AMPLITUDE) & (dfa[AMPLITUDE_COL] < MAX_AMPLITUDE)].copy()
dfa = dfa.loc[dfa[PERIOD_COL].between(PERIOD_MIN, PERIOD_MAX, inclusive="both")].copy()

dfa["month"] = dfa["date"].dt.month

# Keep one row per day
auto_daily = (
    dfa[["date", "month", AMPLITUDE_COL]]
    .drop_duplicates(subset=["date"])
    .sort_values("date")
    .copy()
)

# ======================================================================================
# MERGE DAILY MANUAL + AUTO
# ======================================================================================
df_daily = manual_daily.merge(auto_daily, on=["date", "month"], how="inner").sort_values("date").copy()

# ======================================================================================
# SWEEP THRESHOLDS USING MONTHLY RMSE OF RAW MONTHLY FRACTIONS
# ======================================================================================
thresholds = np.arange(THRESH_START, THRESH_END + 1e-12, THRESH_STEP)

rows = []
auto_frac_by_thr = {}
auto_norm_by_thr = {}
auto_daily_labels_by_thr = {}

for thr in thresholds:
    df_thr = df_daily.copy()
    df_thr["auto_tid_day"] = (df_thr[AMPLITUDE_COL] >= thr).astype(int)

    # DAILY diagnostics for reference only
    y_true = df_thr["manual_tid_day"].values.astype(int)
    y_pred = df_thr["auto_tid_day"].values.astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    f1 = float(f1_score(y_true, y_pred, zero_division=0))
    precision = float(precision_score(y_true, y_pred, zero_division=0))
    recall = float(recall_score(y_true, y_pred, zero_division=0))
    balanced_acc = float(balanced_accuracy_score(y_true, y_pred))
    daily_acc = float(accuracy_score(y_true, y_pred))

    # MONTHLY FRACTIONS
    auto_tid_days = df_thr.groupby("month")["auto_tid_day"].sum().reindex(months, fill_value=0)
    auto_frac = (auto_tid_days / manual_days_in_month).astype(float)

    auto_mean = float(auto_frac.mean())
    auto_norm = auto_frac / auto_mean if auto_mean != 0 else np.full_like(auto_frac.values, np.nan, dtype=float)

    # PRIMARY METRIC: RMSE of RAW monthly fractions
    rmse_raw_monthly_frac = float(np.sqrt(np.mean((manual_frac.values - auto_frac.values) ** 2)))

    # Reference seasonal diagnostics
    seasonal_r_norm = float(np.corrcoef(manual_norm.values, auto_norm)[0, 1]) if np.all(np.isfinite(auto_norm)) else np.nan
    rmse_norm_12mo = float(np.sqrt(np.mean((manual_norm.values - auto_norm) ** 2))) if np.all(np.isfinite(auto_norm)) else np.nan
    seasonal_r_raw = float(np.corrcoef(manual_frac.values, auto_frac.values)[0, 1])

    rows.append({
        "threshold_km": float(thr),

        # PRIMARY SELECTION METRIC
        "rmse_raw_monthly_frac": rmse_raw_monthly_frac,

        # Monthly reference diagnostics
        "seasonal_r_raw_12mo": seasonal_r_raw,
        "seasonal_r_norm_12mo": seasonal_r_norm,
        "rmse_norm_12mo": rmse_norm_12mo,

        # Daily reference diagnostics
        "f1_daily": f1,
        "precision_daily": precision,
        "recall_daily": recall,
        "balanced_accuracy_daily": balanced_acc,
        "accuracy_daily": daily_acc,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),

        # Counts / rates
        "manual_event_rate_daily": float(y_true.mean()),
        "auto_event_rate_daily": float(y_pred.mean()),
        "mean_manual_monthly_frac": manual_mean,
        "mean_auto_monthly_frac": auto_mean,
        "tid_days_manual": int(manual_tid_days.sum()),
        "tid_days_auto": int(auto_tid_days.sum()),
    })

    auto_frac_by_thr[float(thr)] = auto_frac
    auto_norm_by_thr[float(thr)] = auto_norm
    auto_daily_labels_by_thr[float(thr)] = df_thr[["date", "month", "auto_tid_day"]].copy()

sweep_df = pd.DataFrame(rows)

# ======================================================================================
# CHOOSE BEST THRESHOLD
# Rule:
#   1) minimum monthly RMSE(raw fractions)
#   2) if tie, maximum seasonal_r_raw_12mo
#   3) if tie, maximum f1_daily
#   4) if tie, minimum threshold
# ======================================================================================
sweep_sorted = sweep_df.sort_values(
    by=["rmse_raw_monthly_frac", "seasonal_r_raw_12mo", "f1_daily", "threshold_km"],
    ascending=[True, False, False, True]
).reset_index(drop=True)

best_thr = float(sweep_sorted.loc[0, "threshold_km"])
best_rmse_raw = float(sweep_sorted.loc[0, "rmse_raw_monthly_frac"])
best_seasonal_r_raw = float(sweep_sorted.loc[0, "seasonal_r_raw_12mo"])
best_seasonal_r_norm = float(sweep_sorted.loc[0, "seasonal_r_norm_12mo"])
best_rmse_norm = float(sweep_sorted.loc[0, "rmse_norm_12mo"])
best_f1 = float(sweep_sorted.loc[0, "f1_daily"])
best_precision = float(sweep_sorted.loc[0, "precision_daily"])
best_recall = float(sweep_sorted.loc[0, "recall_daily"])
best_bal_acc = float(sweep_sorted.loc[0, "balanced_accuracy_daily"])
best_acc = float(sweep_sorted.loc[0, "accuracy_daily"])

print("\nSweep complete.")
print(
    f"Best threshold = {best_thr:.2f} km | "
    f"RMSE(raw monthly frac) = {best_rmse_raw:.4f} | "
    f"Seasonal r raw = {best_seasonal_r_raw:.3f} | "
    f"Seasonal r norm = {best_seasonal_r_norm:.3f} | "
    f"F1 daily = {best_f1:.3f}"
)

# Save sweep table
csv_out = os.path.join(OUTPUT_DIR, f"{YEAR}_manual_vs_auto_threshold_sweep_monthly_rmse_raw.csv")
print(f"Saving sweep table: {csv_out}")
sweep_df.to_csv(csv_out, index=False)

# ======================================================================================
# BEST THRESHOLD SERIES
# ======================================================================================
best_row = sweep_sorted.loc[0]
auto_frac_best = auto_frac_by_thr[best_thr]
auto_mean_best = float(best_row["mean_auto_monthly_frac"])

# ======================================================================================
# SINGLE FIGURE: 2 PANELS ONLY
#   a) threshold sweep
#   b) raw monthly fractions at best threshold
# ======================================================================================
fig, axes = plt.subplots(2, 1, figsize=FIGSIZE_2PANEL)

# ------------------------------------------------------------------
# Panel a) PRIMARY METRIC vs THRESHOLD
# ------------------------------------------------------------------
ax = axes[0]
ax.plot(
    sweep_df["threshold_km"],
    sweep_df["rmse_raw_monthly_frac"],
    marker="o",
    linewidth=2,
    label="RMSE of raw monthly fractions"
)

ax.axvline(best_thr, linestyle="--", linewidth=2)

ax.set_xlabel("Amplitude Threshold [km]", fontsize=LABEL_FONTSIZE, fontweight="bold")
ax.set_ylabel("RMSE", fontsize=LABEL_FONTSIZE, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)

ax.set_title(
    "a) Threshold Sweep: RMSE of Raw Monthly Fractions",
    fontsize=TITLE_FONTSIZE, fontweight="bold", loc="left", pad=10
)

ax.text(
    0.97, 0.97,
    f"Best thr = {best_thr:.2f} km\n"
    f"Best RMSE = {best_rmse_raw:.4f}",
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=BOX_FONTSIZE,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
)

ax.legend(fontsize=LEGEND_FONTSIZE, loc="best")

# ------------------------------------------------------------------
# Panel b) Raw monthly fractions at selected threshold
# ------------------------------------------------------------------
ax = axes[1]
w = 0.38
ax.bar(
    x - w/2,
    manual_frac.values,
    width=w,
    color="red",
    alpha=BAR_ALPHA,
    edgecolor=EDGE_COLOR,
    label="Manual (2017)"
)
ax.bar(
    x + w/2,
    auto_frac_best.values,
    width=w,
    color="gray",
    alpha=BAR_ALPHA,
    edgecolor=EDGE_COLOR,
    label=f"Automated (2017), thr={best_thr:.2f} km"
)

ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=TICK_FONTSIZE)
ax.set_ylabel("Monthly LSTID-day Fraction", fontsize=LABEL_FONTSIZE, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
ax.set_ylim(0, 0.7)

ax.set_title(
    "b) 2017 Monthly LSTID-day Fraction (Raw)",
    fontsize=TITLE_FONTSIZE, fontweight="bold", loc="left", pad=10
)

ax.text(
    0.97, 0.97,
    f"mean_manual = {manual_mean:.2f}\n"
    f"mean_auto = {auto_mean_best:.2f}",
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=BOX_FONTSIZE,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
)

ax.legend(fontsize=LEGEND_FONTSIZE, loc="lower right")

fig.suptitle(
    f"Manual vs Automated Seasonal Validation ({YEAR})\n"
    f"Manual LSTID day: hours >= {MANUAL_BINARY_HOUR_THRESHOLD:g} hr; "
    f"Automated LSTID day: amp >= {best_thr:.2f} km, {PERIOD_MIN:g}-{PERIOD_MAX:g} hr, R2 >= {R2_MIN:g}",
    fontsize=SUPTITLE_FONTSIZE, fontweight="bold", y=0.98
)

plt.tight_layout(rect=[0, 0, 1, 0.90])

png_out = os.path.join(OUTPUT_DIR, f"{YEAR}_manual_vs_auto_validation_2panel_monthly_rmse_raw.png")
print(f"Saving: {png_out}")
fig.savefig(os.path.join(OUTPUT_DIR, "fig_rmse_thresh.png"), dpi=300, bbox_inches="tight")

if SHOW_FIGURE:
    pass  # (notebook displayed the figure here)
# ======================================================================================
# OPTIONAL: PRINT TOP 20 THRESHOLDS
# ======================================================================================
print("\nTop 20 thresholds by monthly RMSE(raw fractions):")
print(
    sweep_df.sort_values(
        by=["rmse_raw_monthly_frac", "seasonal_r_raw_12mo", "f1_daily", "threshold_km"],
        ascending=[True, False, False, True]
    )[[
        "threshold_km",
        "rmse_raw_monthly_frac",
        "seasonal_r_raw_12mo",
        "seasonal_r_norm_12mo",
        "f1_daily",
        "precision_daily",
        "recall_daily",
        "tid_days_auto"
    ]].head(20).to_string(index=False)
)
plt.close(fig)
print("Done!")
