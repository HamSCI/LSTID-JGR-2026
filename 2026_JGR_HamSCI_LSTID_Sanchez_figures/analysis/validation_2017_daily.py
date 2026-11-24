"""
Day-by-day agreement of the automated and manual 2017 classifications.

Paper Section 4.2.2 (Figure 7). Figure 7 compares monthly fractions; this script prints the
day-level numbers quoted with it:

  - June-August LSTID-active days, automated (33 km threshold) vs manual.
  - December-January counts, automated vs manual.
  - Precision, recall and balanced accuracy on individual days.
  - For summer days flagged only by manual inspection: whether a quality-controlled fit exists,
    and its amplitude compared with the non-LSTID days of the same months.

Manual LSTID day: >= 1 hour of LSTID activity identified by inspection (column `tid` of
2017_data_full_V2.csv), as in Figure 7. Day-level scores use the days present in both
classifications (automated days must pass quality control), as Figure 7 does.

RUN
  python analysis/validation_2017_daily.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import precision_score, recall_score, balanced_accuracy_score
import common as c

man = pd.read_csv(os.path.join(c.DATA_DIR, "2017_data_full_V2.csv"))
man["date"] = pd.to_datetime(man["date"]).dt.normalize()
man = man[man["date"].dt.year == 2017].drop_duplicates("date").set_index("date")
man["manual"] = pd.to_numeric(man["tid"]) >= 1.0

raw = c.load_fits_raw().loc["2017"]
fits = c.load_fits().loc["2017"]
both = man[["manual"]].join(fits[["amplitude_km", "lstid"]], how="inner")
print(f"2017: manual days {len(man)}, quality-controlled automated days {len(fits)}, in both {len(both)}")

summer, midwinter = (6, 7, 8), (12, 1)
for name, months in [("June-August", summer), ("December-January", midwinter)]:
    m_all = man[man.index.month.isin(months)]["manual"].sum()
    m_both = both[both.index.month.isin(months)]["manual"].sum()
    a = both[both.index.month.isin(months)]["lstid"].sum()
    print(f"{name:17s} LSTID-active days: automated {a}, manual {m_all} "
          f"(manual on days that also have a valid fit: {m_both})")

y, p = both["manual"].to_numpy(), both["lstid"].to_numpy()
print(f"\nDaily agreement at 33 km: precision {precision_score(y, p):.2f}, recall {recall_score(y, p):.2f}, "
      f"balanced accuracy {balanced_accuracy_score(y, p):.2f}")

print("\nSummer days flagged only by manual inspection:")
manual_summer = man[man.index.month.isin(summer) & man["manual"]].index
only_manual = [d for d in manual_summer if not (d in both.index and both.loc[d, "lstid"])]
has_fit = [d for d in only_manual if d in fits.index]
no_qc = [d for d in only_manual if d not in fits.index and d in raw.index and pd.notna(raw.loc[d, "amplitude_km"])]
no_fit = [d for d in only_manual if d not in fits.index and d not in no_qc]
print(f"   {len(only_manual)} days: {len(has_fit)} with a quality-controlled fit, "
      f"{len(no_qc)} with a fit that fails quality control, {len(no_fit)} with no fit")
amp_only = fits.loc[has_fit, "amplitude_km"]
non = fits[fits.index.month.isin(summer) & ~fits["lstid"]]["amplitude_km"]
print(f"   amplitude of those with a valid fit: median {amp_only.median():.1f} km "
      f"(IQR {amp_only.quantile(.25):.1f}-{amp_only.quantile(.75):.1f}); "
      f"summer 2017 non-LSTID days: median {non.median():.1f} km "
      f"(IQR {non.quantile(.25):.1f}-{non.quantile(.75):.1f})")
