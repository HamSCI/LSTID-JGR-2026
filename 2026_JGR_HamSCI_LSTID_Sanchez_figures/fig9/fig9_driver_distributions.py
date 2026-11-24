"""
Paper Figure 9: Conditional Distributions of Candidate Drivers on LSTID-Active and Non-LSTID Days

Sanchez et al. (2026), JGR Space Physics, automated LSTID detection with amateur radio.

WHAT IT SHOWS
Estimated probability density functions of (a) AE daily mean, (b) SYM-H daily mean and
(c) polar vortex index, for LSTID-active (red) and non-LSTID (blue) days, 2015-2022.
Distributions are pooled over the year, so separation in (c) largely reflects the shared seasonal
cycle (see paper Section 4.2.4). AE is unavailable for 2020.

Quality control (paper Section 3.5): daily fits with R^2 >= 0.35 and fitted period 1-4.5 hr,
2015-01-01 to 2022-12-31.

INPUTS (in $LSTID_DATA_DIR, default ../data/)
  20140101-20230901_merged.csv
  merra2_pvIndex.csv
  sym_ae.lst

OUTPUT (written next to this script)
  20150101-20221231_event_probability_cpdfs_vertical.png

RUN
  python fig9/fig9_driver_distributions.py      (from the repository root; any directory works)

Extracted from the analysis notebook final_paper_plot_v11.ipynb, cell 6; logic unchanged
except file paths and headless plotting. Printed statistics match the values in the paper.
"""
import os
import matplotlib
matplotlib.use("Agg")   # headless: write files, no display

# Note: the internal class label "Quiet" below means non-LSTID day.
HERE = os.path.dirname(os.path.abspath(__file__))
# Input data folder: set LSTID_DATA_DIR to where the data lives (see data/README.md);
# defaults to ../data/ next to this script.
DATA_DIR = os.environ.get("LSTID_DATA_DIR", os.path.join(HERE, "..", "data"))

# ======================================================================================
# COMBINED FIGURE: LSTID Event Probability Modulation (DECOUPLED DRIVERS)
# 3 panels (3 rows x 1 column) - vertical layout, conditional PDFs only:
# Row 1 (a): AE Index
# Row 2 (b): SYM/H Index
# Row 3 (c): PV Index
# ======================================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# SETTINGS
# =========================
CSV_FILE = os.path.join(DATA_DIR, "20140101-20230901_merged.csv")
SYM_AE_FILE = os.path.join(DATA_DIR, "sym_ae.lst")
PV_FILE = os.path.join(DATA_DIR, "merra2_pvIndex.csv")
OUTPUT_DIR = HERE
START_DATE = "2015-01-01"
END_DATE = "2022-12-31"
R2_COL = "r2"
R2_MIN = 0.35
PERIOD_MIN = 1.0   # paper Sec 3.5: fitted period must lie in 1-4.5 hr
PERIOD_MAX = 4.5
AMPLITUDE_COL = "amplitude_km"
TID_THRESHOLD = 33.0
PROB_BINS = 25
SHOW_FIGURE = True

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================================================================================
# LOAD BASE AMPLITUDE + APPLY GLOBAL QC ONCE
# ======================================================================================
print("Loading base amplitude + applying QC...")
df_base = pd.read_csv(CSV_FILE, parse_dates=True, index_col=0).sort_index()
date_min = pd.to_datetime(START_DATE)
date_max = pd.to_datetime(END_DATE)
df_base = df_base.loc[(df_base.index >= date_min) & (df_base.index <= date_max)].copy()
df_base[R2_COL] = pd.to_numeric(df_base[R2_COL], errors="raise")
pre_n = len(df_base)
df_base = df_base.loc[df_base[R2_COL].ge(R2_MIN) & df_base[R2_COL].notna()].copy()
post_n = len(df_base)
print(f"R2 filter (>={R2_MIN}): kept {post_n}/{pre_n} days")
df_base = df_base.loc[pd.to_numeric(df_base["T_hr"]).between(PERIOD_MIN, PERIOD_MAX, inclusive="both")].copy()
print(f"Period filter ({PERIOD_MIN:g}-{PERIOD_MAX:g} hr): kept {len(df_base)} days")
df_base[AMPLITUDE_COL] = pd.to_numeric(df_base[AMPLITUDE_COL], errors="coerce")
df_base = df_base.dropna(subset=[AMPLITUDE_COL]).copy()
df_base["event_class"] = pd.Series(index=df_base.index, dtype="object")
df_base.loc[df_base[AMPLITUDE_COL] < TID_THRESHOLD, "event_class"] = "Quiet"
df_base.loc[df_base[AMPLITUDE_COL] >= TID_THRESHOLD, "event_class"] = "LSTID"
n_tid_total = int((df_base["event_class"] == "LSTID").sum())
n_quiet_total = int((df_base["event_class"] == "Quiet").sum())
print(f"\nBase after QC: LSTID={n_tid_total}, Quiet={n_quiet_total}, Total={len(df_base)}")

# ======================================================================================
# LOAD DRIVERS
# ======================================================================================
print("\nLoading drivers...")
df_sym_ae = pd.read_csv(
    SYM_AE_FILE,
    sep=r"\s+",
    comment="#",
    header=None,
    names=["year", "day", "hour", "minute", "ae_index", "sym_h"]
)
df_sym_ae["datetime"] = (
    pd.to_datetime(df_sym_ae["year"].astype(str) + df_sym_ae["day"].astype(str), format="%Y%j")
    + pd.to_timedelta(df_sym_ae["hour"], unit="h")
    + pd.to_timedelta(df_sym_ae["minute"], unit="m")
)
df_sym_ae = df_sym_ae.set_index("datetime").sort_index()
df_sym_ae["ae_index"] = df_sym_ae["ae_index"].replace(99999, np.nan)
df_sym_ae["sym_h"] = df_sym_ae["sym_h"].replace(99999, np.nan)
df_sym_ae_daily = df_sym_ae[["ae_index", "sym_h"]].resample("D").mean()
df_pv = pd.read_csv(PV_FILE, header=None, names=["date", "pv_index"])
df_pv["date"] = pd.to_datetime(df_pv["date"], format="%m/%d/%Y", errors="coerce")
df_pv = df_pv.dropna(subset=["date"]).set_index("date").sort_index()

# ======================================================================================
# PROCESS AE INDEX
# ======================================================================================
print("\nProcessing AE Index...")
df_ae = df_base.merge(df_sym_ae_daily[["ae_index"]], left_index=True, right_index=True, how="left")
df_ae = df_ae.dropna(subset=["ae_index"]).copy()
ae_quiet = df_ae.loc[df_ae["event_class"] == "Quiet", "ae_index"].to_numpy()
ae_tid = df_ae.loc[df_ae["event_class"] == "LSTID", "ae_index"].to_numpy()
print(f"AE: Quiet={len(ae_quiet)}, LSTID={len(ae_tid)}")
ae_all = np.concatenate([ae_quiet, ae_tid])
ae_edges = np.linspace(ae_all.min(), ae_all.max(), PROB_BINS + 1)
ae_centers = 0.5 * (ae_edges[:-1] + ae_edges[1:])
ae_pdf_quiet = np.histogram(ae_quiet, bins=ae_edges, density=True)[0] if len(ae_quiet) else np.zeros(PROB_BINS, dtype=float)
ae_pdf_tid = np.histogram(ae_tid, bins=ae_edges, density=True)[0] if len(ae_tid) else np.zeros(PROB_BINS, dtype=float)

# ======================================================================================
# PROCESS SYM-H INDEX
# ======================================================================================
print("Processing SYM-H Index...")
df_symh = df_base.merge(df_sym_ae_daily[["sym_h"]], left_index=True, right_index=True, how="left")
df_symh = df_symh.dropna(subset=["sym_h"]).copy()
symh_quiet = df_symh.loc[df_symh["event_class"] == "Quiet", "sym_h"].to_numpy()
symh_tid = df_symh.loc[df_symh["event_class"] == "LSTID", "sym_h"].to_numpy()
print(f"SYM-H: Quiet={len(symh_quiet)}, LSTID={len(symh_tid)}")
symh_all = np.concatenate([symh_quiet, symh_tid])
symh_edges = np.linspace(symh_all.min(), symh_all.max(), PROB_BINS + 1)
symh_centers = 0.5 * (symh_edges[:-1] + symh_edges[1:])
symh_pdf_quiet = np.histogram(symh_quiet, bins=symh_edges, density=True)[0] if len(symh_quiet) else np.zeros(PROB_BINS, dtype=float)
symh_pdf_tid = np.histogram(symh_tid, bins=symh_edges, density=True)[0] if len(symh_tid) else np.zeros(PROB_BINS, dtype=float)

# ======================================================================================
# PROCESS PV INDEX
# ======================================================================================
print("Processing PV Index...")
df_pv_merge = df_base.merge(df_pv[["pv_index"]], left_index=True, right_index=True, how="left")
df_pv_merge = df_pv_merge.dropna(subset=["pv_index"]).copy()
pv_quiet = df_pv_merge.loc[df_pv_merge["event_class"] == "Quiet", "pv_index"].to_numpy()
pv_tid = df_pv_merge.loc[df_pv_merge["event_class"] == "LSTID", "pv_index"].to_numpy()
print(f"PV: Quiet={len(pv_quiet)}, LSTID={len(pv_tid)}")
pv_all = np.concatenate([pv_quiet, pv_tid])
pv_edges = np.linspace(pv_all.min(), pv_all.max(), PROB_BINS + 1)
pv_centers = 0.5 * (pv_edges[:-1] + pv_edges[1:])
pv_pdf_quiet = np.histogram(pv_quiet, bins=pv_edges, density=True)[0] if len(pv_quiet) else np.zeros(PROB_BINS, dtype=float)
pv_pdf_tid = np.histogram(pv_tid, bins=pv_edges, density=True)[0] if len(pv_tid) else np.zeros(PROB_BINS, dtype=float)

# ======================================================================================
# CREATE 3-PANEL FIGURE (vertical, conditional PDFs only)
# ======================================================================================
print("\nCreating combined figure...")
fig = plt.figure(figsize=(10, 20))
gs = fig.add_gridspec(3, 1, hspace=0.55)

labels = ["a)", "b)", "c)"]
rows = [
    ("AE Index",           "[nT]", ae_centers,   ae_pdf_tid,   ae_pdf_quiet,   ae_tid,   ae_quiet),
    ("SYM/H Index",        "[nT]", symh_centers, symh_pdf_tid, symh_pdf_quiet, symh_tid, symh_quiet),
    ("Polar Vortex Index", "",     pv_centers,   pv_pdf_tid,   pv_pdf_quiet,   pv_tid,   pv_quiet),
]

for r, (name, units, centers, pdf_tid, pdf_quiet, tid, quiet) in enumerate(rows):
    ax = fig.add_subplot(gs[r, 0])
    ax.plot(centers, pdf_tid,   color="red",  lw=1.8, label=f"p({name.split()[0]} | LSTID) (N={len(tid)})")
    ax.plot(centers, pdf_quiet, color="blue", lw=1.8, label=f"p({name.split()[0]} | non-LSTID) (N={len(quiet)})")
    ax.set_xlabel(f"{name} {units}", fontsize=18, fontweight="bold")
    ax.set_ylabel("Probability Density", fontsize=18, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(fontsize=13, loc="upper left" if r == 1 else "upper right")
    ymax = float(np.nanmax([np.nanmax(pdf_tid) if len(pdf_tid) else 0.0, np.nanmax(pdf_quiet) if len(pdf_quiet) else 0.0]))
    ax.set_ylim(0, ymax * 1.10 if ymax > 0 else 1)
    ax.set_title(
        f"{labels[r]} {name}\nEstimated PDFs: LSTID-Active vs Non-LSTID",
        fontsize=20, fontweight="bold", loc="left", pad=10
    )

# MAIN TITLE
fig.suptitle(
    "Estimated Probability Density Functions:\nLSTID-Active vs Non-LSTID Days",
    fontsize=22, fontweight="bold", y=0.99
)
fig.text(
    0.5, 0.95,
    f"{date_min.strftime('%Y-%m-%d')} to {date_max.strftime('%Y-%m-%d')}",
    fontsize=16, fontweight="normal", ha="center", va="top"
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

sDate_str = date_min.strftime("%Y%m%d")
eDate_str = date_max.strftime("%Y%m%d")
png_fpath = os.path.join(OUTPUT_DIR, f"{sDate_str}-{eDate_str}_event_probability_cpdfs_vertical.png")
pdf_fpath = os.path.join(OUTPUT_DIR, f"{sDate_str}-{eDate_str}_event_probability_cpdfs_vertical.pdf")

print(f"\nSaving: {png_fpath}")
fig.savefig(png_fpath, dpi=300, bbox_inches="tight")
print(f"Saving: {pdf_fpath}")
fig.savefig(pdf_fpath, bbox_inches="tight")

if SHOW_FIGURE:
    pass  # (notebook displayed the figure here)
print("\nDone!")
