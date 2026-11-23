# Analysis scripts

Scripts that reproduce the numbers quoted in the paper's text that are not printed by a figure
script. They read the same input data as the figures (`$LSTID_DATA_DIR`, see `../README.md`) and
apply the same quality control (paper Section 3.5). Each prints its results; none writes files.

```bash
export LSTID_DATA_DIR=/path/to/data
python analysis/vortex_seasonal_tests.py    # any script, from the repository folder
```

| Script | Paper section | What it prints |
|---|---|---|
| `vortex_seasonal_tests.py` | 4.2.1, 4.2.4 (Figs 5, 9c) | Amplitude vs polar vortex index after 60- and 180-day means, with and without the day-of-year mean removed; scrambled-year test at 60 and 180 days, with and without 2015; the same test for noon cos(SZA); vortex AUC pooled and within ±15-day calendar windows, 25–35 km thresholds |
| `day_to_day_drivers.py` | 4.2.5 | Season-removed correlation of amplitude with AE, SYM-H, vortex index, F10.7 and skip distance, daily and at 31 days; season-controlled AUC by season |
| `summer_sporadic_e.py` | 4.2.5 | May–August high vs low sporadic-E days: skip distance, amplitude, valid-fit fraction, LSTID fraction, permutation p values; amplitude as a fraction of skip distance by month |
| `validation_2017_daily.py` | 4.2.2 (Fig 7) | 2017 summer and midwinter counts, automated vs manual; daily precision, recall, balanced accuracy; what the manual-only summer days look like |
| `qc_and_skip_windows.py` | 3.5, 4.3.1 (Figs 8, 10) | Days kept by each quality-control step; skip distance vs F10.7 at 180 and 365 days; annual-cycle gain and independent samples per window |

`common.py` holds the loaders, the day-of-year climatology (±15 days, as in Figure 10e), the
moving-block bootstrap (as in Figures 5 and 10) and the AUC helpers the scripts share.

Runtimes are under a minute each, except `vortex_seasonal_tests.py` (14,833 derangements and
the bootstrap, a few minutes).
