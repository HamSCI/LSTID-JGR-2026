"""
Do the candidate drivers explain day-to-day changes in fitted amplitude?

Paper Section 4.2.5: "once the seasonal cycle was removed, geomagnetic activity showed a weak but
detectable association, whereas the polar vortex index showed none." This script prints:

  1. Correlation of amplitude with each driver after removing the day-of-year mean (+-15 d) from
     both series. Rows are daily values and 31-day running means; the 31-day rows are repeated
     with the year-to-year (solar-cycle) variation also removed, by subtracting a 365-day running
     mean. Comparing the rows shows which timescale a correlation comes from.
  2. Season-controlled AUC: how well each driver separates LSTID-active from non-LSTID days when
     days are compared only within +-15-day calendar windows, overall and by season.

Drivers: AE averaged over each day's fit window, AE daily mean (as in Figure 9), AE daily
maximum, SYM-H daily mean, polar vortex index, F10.7, and minimum skip distance from the same fit.
AE is unavailable for 2020. Quality control as in Figure 9 (paper Section 3.5).

RUN
  python analysis/day_to_day_drivers.py
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import common as c

fits = c.load_fits()
geo = c.load_ae_symh_5min()
grid = c.calendar()


def fit_window_mean(row):
    """Mean 5-minute AE between the fit's start and end times."""
    try:
        t0 = pd.to_datetime(row["fitStart"], format="%m/%d/%Y %H:%M")
        t1 = pd.to_datetime(row["fitEnd"], format="%m/%d/%Y %H:%M")
    except (TypeError, ValueError):
        return np.nan
    return geo["ae"].loc[t0:t1].mean()


d = pd.DataFrame(index=grid)
d["amp"] = fits["amplitude_km"]
d["AE, fit window"] = fits.apply(fit_window_mean, axis=1)
d["AE, daily mean"] = geo["ae"].resample("D").mean()
d["AE, daily max"] = geo["ae"].resample("D").max()
d["SYM-H, daily mean"] = geo["symh"].resample("D").mean()
d["polar vortex index"] = c.load_pv()
d["F10.7"] = c.load_f107()
d["min skip distance"] = fits["min_combined_fit"]
drivers = [k for k in d.columns if k != "amp"]
anom = d.apply(c.deseason)


def corr(x, y):
    ok = x.notna() & y.notna()
    return pearsonr(x[ok], y[ok])[0], spearmanr(x[ok], y[ok])[0], int(ok.sum())


def remove_interannual(s):
    return s - s.rolling(365, center=True, min_periods=60).mean()


print("1. Amplitude vs driver, day-of-year mean removed from both (Pearson r / Spearman rho / N)")
print(f"   {'driver':20s} {'daily':>22s} {'31-day means':>22s} {'31-day, no year-to-year':>25s}")
for k in drivers:
    a, x = anom["amp"], anom[k]
    daily = corr(a, x)
    sm = corr(c.running_mean(a, 31).where(a.notna()), c.running_mean(x, 31).where(a.notna()))
    a2, x2 = remove_interannual(a), remove_interannual(x)
    sm2 = corr(c.running_mean(a2, 31).where(a.notna()), c.running_mean(x2, 31).where(a.notna()))
    print(f"   {k:20s} " + "  ".join(f"{p:+.2f} / {s:+.2f} / {n:4d}" for p, s, n in (daily, sm, sm2)))

print("\n   Annual means (year-to-year variation behind the 31-day rows):")
yr = d.groupby(d.index.year)[["amp", "AE, daily mean", "F10.7", "min skip distance"]].mean()
print(yr.round(1).to_string().replace("\n", "\n   "))

print("\n2. Season-controlled AUC, LSTID-active vs non-LSTID days (0.5 = no separation;")
print("   below 0.5 means LSTID-active days have lower values)")
seasons = {"all": range(1, 13), "DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
print(f"   {'driver':20s} " + " ".join(f"{s:>6s}" for s in seasons))
active = fits["lstid"].reindex(grid)
for k in ["AE, fit window", "AE, daily mean", "AE, daily max", "SYM-H, daily mean", "polar vortex index"]:
    vals = []
    for months in seasons.values():
        m = active.notna() & d[k].notna() & d.index.month.isin(months)
        vals.append(c.auc_within_calendar(d.loc[m, k], active[m].astype(bool)))
    print(f"   {k:20s} " + " ".join(f"{v:6.3f}" for v in vals))
