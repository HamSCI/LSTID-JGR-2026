"""
Quality-control counts, and the choice of smoothing window for skip distance vs F10.7.

Paper Sections 3.5 and 4.3.1 (Figures 8 and 10).

  1. How many days each quality-control step keeps (paper Section 3.5), and the LSTID-active /
     non-LSTID split at 33 km (Figure 8).
  2. Minimum skip distance vs F10.7 after 180-day and 365-day running means, each with a 95%
     moving-block bootstrap interval (block = window), built as in Figure 10.
  3. How much of an annual cycle survives each running mean (the gain of a boxcar of that
     length at a period of 365.25 days), and the number of independent samples each window
     leaves over the record.

RUN
  python analysis/qc_and_skip_windows.py
"""
import numpy as np
import common as c

# ---------------------------------------------------------------- 1. QC counts
raw = c.load_fits_raw()
has_fit = raw["r2"].notna()
r2_ok = raw["r2"].ge(c.R2_MIN)
per_ok = raw["T_hr"].between(c.PERIOD_MIN, c.PERIOD_MAX, inclusive="both")
print("1. Quality control, 2015-01-01 to 2022-12-31")
print(f"   calendar days {len(c.calendar())}; days with a fit {has_fit.sum()}")
print(f"   R^2 >= {c.R2_MIN}: {r2_ok.sum()}")
print(f"   ... and period {c.PERIOD_MIN:g}-{c.PERIOD_MAX:g} hr: {(r2_ok & per_ok).sum()} "
      f"(period limit removes {(r2_ok & ~per_ok).sum()})")
fits = c.load_fits()
print(f"   LSTID-active (>= {c.TID_THRESHOLD:g} km) {fits['lstid'].sum()}, non-LSTID {(~fits['lstid']).sum()}")

# ---------------------------------------------------------------- 2. skip vs F10.7
frame = c.fit_day_frame(fits, f107=c.load_f107())
ms, fx = frame["min_combined_fit"], frame["f107"]
pair_days = ms.notna() & fx.notna()
print("\n2. Minimum skip distance vs F10.7 (95% block-bootstrap CI, block = window)")
for w in (180, 365):
    r, lo, hi, n = c.block_bootstrap_correlation(c.running_mean(ms, w).where(pair_days),
                                                 c.running_mean(fx, w).where(pair_days), block_length=w)
    print(f"   {w}-day mean: r = {r:+.3f}  CI [{lo:+.3f}, {hi:+.3f}]  N = {n}")

# ---------------------------------------------------------------- 3. window properties
print("\n3. Annual cycle passed by a running mean, and independent samples over the record")
span = (frame.index[-1] - frame.index[0]).days + 1
for w in (60, 180, 365):
    gain = abs(np.sin(np.pi * w / 365.25) / (w * np.sin(np.pi / 365.25)))
    print(f"   {w:3d}-day mean: passes {100 * gain:4.1f}% of the annual-cycle amplitude; "
          f"{span / w:4.1f} window lengths in {span} days")
