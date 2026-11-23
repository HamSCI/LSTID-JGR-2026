"""
Does summer sporadic-E suppress the detection of LSTID-like structure?

Paper Section 4.2.5. In summer, sporadic-E can return the 14 MHz signal below the F2 layer, so
the skip-distance edge may be set by a layer that carries no F-region wave structure. The test:
among May-August days, compare days of high and low sporadic-E occurrence after removing the
day-of-year mean from the occurrence series. If sporadic-E hid LSTIDs, high-sporadic-E days
would show smaller amplitudes, fewer valid fits, or fewer LSTID-active days.

  - Sporadic-E occurrence: percent of daytime (13-23 UTC) ionograms with foEs >= 4 MHz across
    six GIRO stations, as in Figure 11. High / low = more / less sporadic-E than the day-of-year
    mean (anomaly above / at or below zero). Every May-August day is in one group.
  - Valid fit: the day passes the quality control of paper Section 3.5.
  - p values: permutation tests that shuffle the high/low labels within each calendar month of
    each year (10000 permutations, two-sided, difference of means).

Also prints the monthly climatology of amplitude as a fraction of minimum skip distance (the
normalized-amplitude check in the same section).

RUN
  python analysis/summer_sporadic_e.py
"""
import numpy as np
import pandas as pd
import common as c

SUMMER = (5, 6, 7, 8)
N_PERM, SEED = 10000, 42

fits = c.load_fits()
es, _ = c.load_giro_daily()
grid = c.calendar()
d = pd.DataFrame(index=grid)
d["es"] = es
d["es_anom"] = c.deseason(d["es"])
d["valid"] = d.index.isin(fits.index).astype(float)
d["amp"] = fits["amplitude_km"]
d["skip"] = fits["min_combined_fit"]
d["lstid"] = fits["lstid"].astype(float)
s = d[d.index.month.isin(SUMMER) & d["es_anom"].notna()].copy()
s["high"] = s["es_anom"] > 0
s["stratum"] = s.index.year * 100 + s.index.month
print(f"May-August days with sporadic-E data: {len(s)} (high {s['high'].sum()}, low {(~s['high']).sum()})")

rng = np.random.default_rng(SEED)
groups = [np.flatnonzero(s["stratum"].to_numpy() == k) for k in np.unique(s["stratum"])]


def diff(values, labels):
    ok = ~np.isnan(values)
    return values[ok & labels].mean() - values[ok & ~labels].mean()


def perm_p(col):
    v, lab = s[col].to_numpy(float), s["high"].to_numpy()
    obs = diff(v, lab)
    count = 0
    for _ in range(N_PERM):
        perm = lab.copy()
        for g in groups:
            perm[g] = rng.permutation(perm[g])
        count += abs(diff(v, perm)) >= abs(obs)
    return obs, (count + 1) / (N_PERM + 1)


print(f"\n{'quantity':34s} {'high Es':>9s} {'low Es':>9s} {'difference':>11s} {'p':>7s}")
for col, label, fmt in [("skip", "minimum skip distance [km]", "{:9.0f}"),
                        ("amp", "fitted amplitude [km]", "{:9.1f}"),
                        ("valid", "fraction of days with a valid fit", "{:9.2f}"),
                        ("lstid", "fraction LSTID-active (of fits)", "{:9.2f}")]:
    hi, lo = s.loc[s["high"], col].mean(), s.loc[~s["high"], col].mean()
    obs, p = perm_p(col)
    print(f"{label:34s} {fmt.format(hi)} {fmt.format(lo)} {obs:+11.2f} {p:7.3f}")

print("\nAmplitude as a fraction of minimum skip distance, monthly mean over 2015-2022 (x 100):")
ratio = (fits["amplitude_km"] / fits["min_combined_fit"]) * 100
clim = ratio.groupby(ratio.index.month).mean()
print("   " + "  ".join(f"{m:>4d}" for m in range(1, 13)))
print("   " + "  ".join(f"{clim[m]:4.2f}" for m in range(1, 13)))
