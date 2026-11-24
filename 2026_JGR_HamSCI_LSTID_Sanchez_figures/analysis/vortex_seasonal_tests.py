"""
Is the amplitude-vortex relationship more than a shared seasonal cycle?

Paper Sections 4.2.1 and 4.2.4 (Figures 5 and 9c). Both fitted amplitude and the polar vortex
index peak in winter, so statistics pooled over the year partly measure agreement of the two
annual cycles. This script prints the tests the text quotes:

  1. Amplitude vs vortex after a 60-day running mean (Figure 5d), and the same with the
     day-of-year mean removed from each series before smoothing, at 60 and 180 days, with 95%
     moving-block bootstrap intervals (block = window).
  2. Scrambled-year test: each year's smoothed amplitude paired with a different year's smoothed
     vortex index, over every derangement of the years. Reports the median r and the fraction of
     derangements with r at least the observed value, at 60 and 180 days, with all years and
     without 2015 (whose winter is only January-February, has the weakest vortex and the highest
     F10.7 of the record). Supports "covariation beyond the seasonal cycle could not be
     established" (Section 4.2.1): agreement between years at 180 days rests on that one winter.
  3. Noon cos(solar zenith angle), a pure function of the date, on the same 60-day test, with and
     without the day-of-year mean removed.
  4. Separation of LSTID-active and non-LSTID days by the vortex index (Figure 9c) as an AUC,
     pooled over the year and within +-15-day calendar windows, at thresholds 25-35 km.

Quality control and series construction as in Figure 5 (paper Section 3.5): the vortex index is
kept on days with a quality-controlled fit, smoothed values are compared on those days only, and
the record ends 2022-07-30 with the vortex index. The AUCs use all quality-controlled days, as
Figure 9 does.

RUN
  python analysis/vortex_seasonal_tests.py
"""
import itertools
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import common as c

fits = c.load_fits()
frame = c.fit_day_frame(fits, pv=c.load_pv(), cos_sza=c.load_cos_sza())
grid = frame.index
amp, pv, cos_sza = frame["amplitude_km"], frame["pv"], frame["cos_sza"]
pair_days = amp.notna() & pv.notna()           # smoothed values are compared on these days only
print(f"Daily grid {grid[0].date()} to {grid[-1].date()}: {pair_days.sum()} days with amplitude and vortex")


def on_pair_days(s):
    return s.where(pair_days)

# ---------------------------------------------------------------- 1. smoothed correlations
print("\n1. Amplitude vs polar vortex index (95% block-bootstrap CI, block = window)")
print(f"   daily r = {pearsonr(amp[pair_days], pv[pair_days])[0]:+.3f}")
for label, a, p in [("raw", amp, pv), ("day-of-year mean removed", c.deseason(amp), c.deseason(pv))]:
    for w in (60, 180):
        r, lo, hi, n = c.block_bootstrap_correlation(on_pair_days(c.running_mean(a, w)),
                                                     on_pair_days(c.running_mean(p, w)), block_length=w)
        print(f"   {label:26s} {w:3d}-day mean: r = {r:+.3f}  CI [{lo:+.2f}, {hi:+.2f}]  N = {n}")

# ---------------------------------------------------------------- 2. scrambled years
print("\n2. Scrambled-year test (each year paired with a different year)")
doy = c.doy_index(grid)


def scrambled_year_test(w, drop=()):
    a_sm, p_sm = on_pair_days(c.running_mean(amp, w)), on_pair_days(c.running_mean(pv, w))
    years = [y for y in sorted(set(grid.year)) if y not in drop]
    A = np.full((len(years), 365), np.nan)
    P = np.full((len(years), 365), np.nan)
    for i, y in enumerate(years):
        m = (grid.year == y) & ~((grid.month == 2) & (grid.day == 29))
        A[i, doy[m] - 1] = a_sm[m].to_numpy()
        P[i, doy[m] - 1] = p_sm[m].to_numpy()

    def paired_r(perm):
        x, y = A.ravel(), P[list(perm)].ravel()
        ok = ~np.isnan(x) & ~np.isnan(y)
        return pearsonr(x[ok], y[ok])[0]

    n = len(years)
    r_obs = paired_r(range(n))
    r_scr = np.array([paired_r(p) for p in itertools.permutations(range(n))
                      if all(p[i] != i for i in range(n))])
    label = "all years" if not drop else "without " + ", ".join(map(str, drop))
    print(f"   {w:3d}-day mean, {label:13s}: observed r = {r_obs:+.3f}; {len(r_scr)} derangements: "
          f"median r = {np.median(r_scr):+.3f}, fraction with r >= observed p = {np.mean(r_scr >= r_obs):.3f}")


for w in (60, 180):
    for drop in ((), (2015,)):
        scrambled_year_test(w, drop)

# ---------------------------------------------------------------- 3. cos(SZA)
print("\n3. Amplitude vs noon cos(SZA), 60-day means (block-bootstrap CI)")
for label, a, s in [("raw", amp, cos_sza), ("day-of-year mean removed", c.deseason(amp), c.deseason(cos_sza))]:
    if label != "raw" and np.nanstd(c.deseason(cos_sza)) < 1e-3:
        print(f"   {label:26s}: cos(SZA) has no variation left once its day-of-year mean is removed")
        continue
    r, lo, hi, n = c.block_bootstrap_correlation(on_pair_days(c.running_mean(a, 60)),
                                                 on_pair_days(c.running_mean(s, 60)), block_length=60)
    print(f"   {label:26s}: r = {r:+.3f}  CI [{lo:+.2f}, {hi:+.2f}]")

# ---------------------------------------------------------------- 4. AUC
print("\n4. Vortex index on LSTID-active vs non-LSTID days, as AUC (0.5 = no separation)")
print("   threshold  N active  AUC pooled  AUC within +-15 d calendar window")
fits = c.load_fits()                            # as Figure 9
pv = c.load_pv().reindex(c.calendar())
cos_sza = c.load_cos_sza().reindex(c.calendar())
grid = c.calendar()
for thr in (25.0, 30.0, 33.0, 35.0):
    act = (fits["amplitude_km"] >= thr).reindex(grid)
    ok = act.notna() & pv.notna()
    v, a = pv[ok], act[ok].astype(bool)
    print(f"   {thr:5.0f} km   {a.sum():6d}    {c.auc(v[a], v[~a]):.3f}       {c.auc_within_calendar(v, a):.3f}")
act = fits["lstid"].reindex(grid)
ok = act.notna() & cos_sza.notna()
v, a = cos_sza[ok], act[ok].astype(bool)
print(f"   noon cos(SZA), 33 km: pooled {c.auc(v[a], v[~a]):.3f}, within calendar {c.auc_within_calendar(v, a):.3f}")
