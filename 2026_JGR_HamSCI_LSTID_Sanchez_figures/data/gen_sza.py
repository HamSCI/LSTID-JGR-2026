"""
Generate a daily solar zenith angle (SZA) dataset at solar noon.

Location: geographic center of the contiguous US (39.83 N, 98.58 W).
Method:   for each day, sample pysolar altitude every 10 minutes over the
          full 24 h window (UTC) and record the minimum SZA (= maximum
          solar elevation, i.e. local solar noon).
Output:   $LSTID_DATA_DIR/sza_daily.csv  — columns: date, sza_noon_deg
"""

import datetime
import os
import numpy as np
import pandas as pd
from pysolar.solar import get_altitude

LAT =  39.8333   # deg N — geographic center of contiguous US
LON = -98.5833   # deg W

START = datetime.date(2014, 1, 1)
END   = datetime.date(2023, 9, 1)

SAMPLE_MINUTES = 10   # resolution for finding solar noon within each day

# Output folder: LSTID_DATA_DIR if set (see README.md), otherwise this script's folder.
DATA_DIR = os.environ.get("LSTID_DATA_DIR", os.path.dirname(os.path.abspath(__file__)))
OUT_FILE = os.path.join(DATA_DIR, "sza_daily.csv")

records = []
date = START
while date <= END:
    # Build UTC datetimes for this day at SAMPLE_MINUTES intervals
    altitudes = []
    for minute_offset in range(0, 1440, SAMPLE_MINUTES):
        dt = datetime.datetime(
            date.year, date.month, date.day,
            minute_offset // 60, minute_offset % 60, 0,
            tzinfo=datetime.timezone.utc
        )
        alt = get_altitude(LAT, LON, dt)
        altitudes.append(alt)

    max_alt = max(altitudes)
    sza_noon = 90.0 - max_alt   # SZA at solar noon

    records.append({"date": date.isoformat(), "sza_noon_deg": round(sza_noon, 4)})
    date += datetime.timedelta(days=1)

df = pd.DataFrame(records).set_index("date")
df.index = pd.to_datetime(df.index)
df.to_csv(OUT_FILE)

print(f"Saved {OUT_FILE}  ({len(df)} days, {df.index.min().date()} – {df.index.max().date()})")
print(df.describe())
