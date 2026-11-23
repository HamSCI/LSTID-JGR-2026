"""
Shared loaders and statistics for the analysis scripts in this folder.

The figure scripts are deliberately standalone. The analysis scripts share this module because
they all start from the same daily series (quality-controlled fits, polar vortex index, AE, SYM-H,
F10.7, noon solar zenith angle, GIRO sporadic-E and foF2), and every one of them removes a seasonal
cycle the same way. Each loader repeats the reading and quality-control choices of the figure
script named in its docstring, so a number printed here and the matching figure come from the
same data.

Data folder: $LSTID_DATA_DIR, default ../data/ (see ../README.md and ../data/README.md).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.environ.get("LSTID_DATA_DIR", os.path.join(HERE, "..", "data"))

START_DATE, END_DATE = "2015-01-01", "2022-12-31"
R2_MIN = 0.35                 # paper Section 3.5
PERIOD_MIN, PERIOD_MAX = 1.0, 4.5
MIN_AMPLITUDE, MAX_AMPLITUDE = 0.0, 90.0
MIN_SKIP_MIN, MIN_SKIP_MAX = 0.0, 1600.0
TID_THRESHOLD = 33.0          # LSTID-active day: amplitude >= 33 km (Figure 7)
CLIM_HALF_WINDOW = 15         # +-days for the day-of-year climatology (as in Figure 10e)

GIRO_STATIONS = ["BC840", "AU930", "EG931", "AL945", "WP937", "MHJ45"]
GIRO_HOURS = (13, 23)         # daytime UTC window (Figure 11)
ES_THRESH = 4.0               # foEs >= 4 MHz counts as sporadic-E occurrence (Figure 11)


def fit_day_frame(fits, **series):
    """Daily frame built the way Figures 5 and 10 build it: the drivers are kept only on days with
    a quality-controlled fit (inner merge), then the frame is put on a daily grid so that running
    means skip the missing days."""
    df = fits[["amplitude_km", "min_combined_fit"]].copy()
    for name, s in series.items():
        df = df.join(s.rename(name), how="inner")
    return df.sort_index().asfreq("D")


def calendar():
    """Daily DatetimeIndex over the study interval."""
    return pd.date_range(START_DATE, END_DATE, freq="D")


def load_fits_raw():
    """All daily fit rows in the study interval, before quality control (numeric columns parsed)."""
    df = pd.read_csv(os.path.join(DATA_DIR, "20140101-20230901_merged.csv"),
                     parse_dates=True, index_col=0).sort_index()
    df = df.loc[START_DATE:END_DATE].copy()
    for c in ["r2", "T_hr", "amplitude_km", "min_combined_fit"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_fits(skip_range=False):
    """Quality-controlled daily fits, as in Figures 8 and 9: R^2 >= 0.35, period 1-4.5 hr,
    0 < amplitude < 90 km. With skip_range=True, also 0 < minimum skip distance <= 1600 km, as in
    Figures 5, 10 and 11. Adds `lstid` (amplitude >= 33 km)."""
    df = load_fits_raw()
    if skip_range:
        df = df.loc[(df["min_combined_fit"] > MIN_SKIP_MIN) & (df["min_combined_fit"] <= MIN_SKIP_MAX)]
    df = df.loc[df["r2"].ge(R2_MIN)]
    df = df.loc[df["T_hr"].between(PERIOD_MIN, PERIOD_MAX, inclusive="both")]
    df = df.loc[(df["amplitude_km"] > MIN_AMPLITUDE) & (df["amplitude_km"] < MAX_AMPLITUDE)].copy()
    df["lstid"] = df["amplitude_km"] >= TID_THRESHOLD
    return df


def load_pv():
    """Daily polar vortex index (merra2_pvIndex.csv), as read by Figures 5 and 9."""
    pv = pd.read_csv(os.path.join(DATA_DIR, "merra2_pvIndex.csv"), header=None, names=["date", "pv"])
    pv["date"] = pd.to_datetime(pv["date"], format="%m/%d/%Y", errors="coerce")
    pv = pv.dropna(subset=["date"]).drop_duplicates("date").set_index("date").sort_index()
    return pv["pv"]


def load_ae_symh_5min():
    """5-minute AE and SYM-H (sym_ae.lst), fill values removed, as read by Figure 9."""
    d = pd.read_csv(os.path.join(DATA_DIR, "sym_ae.lst"), sep=r"\s+", comment="#", header=None,
                    names=["year", "day", "hour", "minute", "ae", "symh"])
    d.index = (pd.to_datetime(d["year"].astype(str) + d["day"].astype(str), format="%Y%j")
               + pd.to_timedelta(d["hour"], unit="h") + pd.to_timedelta(d["minute"], unit="m"))
    return d[["ae", "symh"]].replace(99999, np.nan).sort_index()


def load_f107():
    """Daily F10.7 from the OMNI file, as read by Figure 10 (values outside 0-200 set missing)."""
    d = pd.read_csv(os.path.join(DATA_DIR, "omni2_daily_3hWFoJxMID.lst"), sep=r"\s+", header=None,
                    names=["year", "day", "hour", "proton_density", "f107"])
    d.index = pd.to_datetime(d["year"].astype(str) + d["day"].astype(str), format="%Y%j")
    f = pd.to_numeric(d["f107"])
    return f.where((f > 0) & (f <= 200))


def load_cos_sza():
    """Daily cosine of the noon solar zenith angle at 39.83 N, 98.58 W (sza_daily.csv, gen_sza.py)."""
    d = pd.read_csv(os.path.join(DATA_DIR, "sza_daily.csv"), parse_dates=["date"], index_col="date")
    return np.cos(np.deg2rad(d["sza_noon_deg"]))


def load_giro_daily():
    """Daily sporadic-E occurrence (percent of daytime ionograms with foEs >= 4 MHz, missing foEs
    counted as no occurrence) and daily mean foF2, over the six GIRO stations, as in Figure 11."""
    rows = []
    t0, t1 = pd.Timestamp(START_DATE), pd.Timestamp(END_DATE) + pd.Timedelta(days=1)
    for stn in GIRO_STATIONS:
        for yr in range(2015, 2023):
            with open(os.path.join(DATA_DIR, "giro", f"giro_{stn}_{yr}.txt")) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    p = line.split()
                    if len(p) < 24 or float(p[1]) < 0:
                        continue
                    t = pd.Timestamp(p[0]).tz_convert(None)
                    if t.hour < GIRO_HOURS[0] or t.hour >= GIRO_HOURS[1] or t < t0 or t >= t1:
                        continue
                    rows.append((t.normalize(),
                                 np.nan if p[12] == "---" else float(p[12]),
                                 np.nan if p[2] == "---" else float(p[2])))
    g = pd.DataFrame(rows, columns=["date", "foEs", "foF2"])
    es = g.groupby("date")["foEs"].apply(lambda s: (s.fillna(-9) >= ES_THRESH).sum() / len(s)) * 100.0
    f2 = g.groupby("date")["foF2"].mean()
    return es, f2


def doy_index(idx):
    """Day of year 1-365, with 29 February merged into 28 February (as in Figure 10e)."""
    doy = idx.dayofyear.to_numpy().copy()
    doy[idx.is_leap_year & (doy > 59)] -= 1
    return np.clip(doy, 1, 365)


def doy_climatology(series, half_window=CLIM_HALF_WINDOW):
    """Day-of-year mean of a daily series, averaged over all years within +-half_window days
    (circular), returned on the series' own index. Same estimator as Figure 10e."""
    vals = series.to_numpy(dtype=float)
    doy = doy_index(series.index)
    clim = np.full(365, np.nan)
    for k in range(1, 366):
        dist = np.minimum(np.abs(doy - k), 365 - np.abs(doy - k))
        m = (dist <= half_window) & ~np.isnan(vals)
        if m.any():
            clim[k - 1] = vals[m].mean()
    return pd.Series(clim[doy - 1], index=series.index)


def deseason(series):
    """Series minus its day-of-year climatology."""
    return series - doy_climatology(series)


def running_mean(series, window):
    """Centered running mean with min_periods = window // 2, as in the figure scripts."""
    return series.rolling(window, center=True, min_periods=max(1, window // 2)).mean()


def block_bootstrap_correlation(x, y, *, block_length, n_bootstrap=2000, seed=42):
    """Pearson r with a 95% moving-block bootstrap interval; same estimator as Figures 5 and 10.
    Returns (r, ci_low, ci_high, n)."""
    ok = x.notna() & y.notna()
    xa, ya = x[ok].to_numpy(), y[ok].to_numpy()
    n = len(xa)
    r_obs = pearsonr(xa, ya)[0]
    rng = np.random.default_rng(seed)
    starts = np.arange(0, n - block_length + 1)
    boot = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        idx = []
        while len(idx) < n:
            s = rng.choice(starts)
            idx.extend(range(s, s + block_length))
        idx = np.asarray(idx[:n])
        boot[b] = pearsonr(xa[idx], ya[idx])[0]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return r_obs, lo, hi, n


def auc(pos, neg):
    """Probability that a random value from `pos` exceeds one from `neg` (ties count one half)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().to_numpy()
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def auc_within_calendar(values, active, half_window=CLIM_HALF_WINDOW):
    """Season-controlled AUC: each LSTID-active day is compared only with non-LSTID days (any year)
    whose day of year is within +-half_window days. Pooled over all such pairs; ties count one half."""
    v = values.to_numpy(dtype=float)
    a = active.to_numpy(dtype=bool)
    ok = ~np.isnan(v)
    v, a, doy = v[ok], a[ok], doy_index(values.index[ok])
    wins = pairs = 0.0
    for i in np.flatnonzero(a):
        dist = np.minimum(np.abs(doy - doy[i]), 365 - np.abs(doy - doy[i]))
        neg = v[(~a) & (dist <= half_window)]
        wins += (v[i] > neg).sum() + 0.5 * (v[i] == neg).sum()
        pairs += len(neg)
    return wins / pairs if pairs else np.nan
