# Figures and data for Sanchez et al. (2026), JGR Space Physics

Scripts and input data that reproduce the figures in:

> Sanchez, D. F., Frissell, N. A., et al. (2026). *Multi-Year Climatology of Large-Scale Traveling
> Ionospheric Disturbances Observed with High Frequency Amateur Radio Receiving Networks Using a
> Novel Automated Detection Algorithm.* Journal of Geophysical Research: Space Physics. (in preparation)

Each paper figure has its own folder. Run a script and the figure is written next to it.

| Folder | Paper figure | Script |
|---|---|---|
| `fig1/` | Example day of amateur radio observations (15 Dec 2018) | not generated here, see `fig1/README.md` |
| `fig2/` | Pre-processing steps for one day | not generated here, see `fig2/README.md` |
| `fig3/` | Detection and sinusoidal fit for one day | not generated here, see `fig3/README.md` |
| `fig4/` | Daily fit parameters, 2015-2022 | `fig4_fit_parameters.py` |
| `fig5/` | Fitted amplitude vs polar vortex index | `fig5_amplitude_vs_polar_vortex.py` |
| `fig6/` | Monthly climatology of amplitude, polar vortex, AE, SYM-H | `fig6_seasonal_phase.py` |
| `fig7/` | Manual vs automated validation, 2017 | `fig7_manual_validation_2017.py` |
| `fig8/` | Monthly LSTID-active and non-LSTID days | `fig8_monthly_lstid_counts.py` |
| `fig9/` | Driver distributions on LSTID-active vs non-LSTID days | `fig9_driver_distributions.py` |
| `fig10/` | Minimum skip distance vs F10.7, with seasonal residual | `fig10_minskip_vs_f107.py` |
| `fig11/` | Minimum skip distance vs sporadic-E and foF2, by season | `fig11_minskip_vs_ionosonde.py` |

Every script starts with a docstring describing the figure, its inputs, processing, and output.

Numbers quoted in the text that no figure script prints (season-removed correlations, the
scrambled-year and sporadic-E tests, the 2017 day-level scores, quality-control counts) are
reproduced by the scripts in `analysis/`; see `analysis/README.md`.

## Running

```bash
pip install -r requirements.txt
export LSTID_DATA_DIR=/path/to/data    # folder holding the input data (see Data below)
python fig4/fig4_fit_parameters.py     # one figure
./run_all.sh                           # all figures, then the analysis scripts
```

Scripts read their inputs from the folder named by the `LSTID_DATA_DIR` environment variable. If it
is not set, they fall back to `data/` in this repository, so they run from any directory. Each script
also prints the statistics quoted in the paper (correlations, confidence intervals, sample sizes).

Tested with Python 3.13.9, numpy 2.3.5, pandas 2.3.3, matplotlib 3.10.6, scipy 1.16.3. With these
versions the outputs are byte-identical to the published figures.

## Data

The input data are not stored in this repository. Put them in one folder (with the GIRO files in
a `giro/` subfolder) and point `LSTID_DATA_DIR` at it. TODO (author): say where to obtain the data.
See `data/README.md` for file formats and provenance and
`DATA_LICENSES.md` for the terms under which each dataset is redistributed.

## Common processing

All figures use the same quality control on the daily fits (paper Section 3.5): goodness of fit
R^2 >= 0.35 and fitted period between 1 and 4.5 hours, over 2015-01-01 to 2022-12-31. LSTID-active
days are quality-controlled days with fitted amplitude >= 33 km (threshold from Figure 7).

## License

Code: TODO (author to choose, e.g. MIT or BSD-3-Clause). Data: see `DATA_LICENSES.md`.
