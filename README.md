# LSTID-JGR-2026

Code to reproduce the figures of *Multi-Year Climatology of Large-Scale Traveling Ionospheric
Disturbances Observed with High Frequency Amateur Radio Receiving Networks Using a Novel
Automated Detection Algorithm* (Sanchez et al., in preparation for *JGR: Space Physics*).

**Status: being assembled.** The figure scripts are being moved here from the authors' working
code. This README will gain installation and per-figure instructions as they arrive, and the
version cited by the paper will be archived on Zenodo.

## How this fits together

- **LSTID detection**: [`HamSCI/hamsci_LSTID_detection`](https://github.com/HamSCI/hamsci_LSTID_detection),
  released separately (release `V0.2.0_Sanchez_et_al_2026_MS`, Zenodo DOI
  [10.5281/zenodo.20186618](https://doi.org/10.5281/zenodo.20186618)). It turns daily
  HamSCI amateur radio spot files from the [CEDAR Madrigal database](http://cedar.openmadrigal.org)
  into per-day LSTID detections and sinusoidal fit parameters.
- **This repository**: starts from that output plus external geophysical index data and produces
  the paper's figures and the statistics quoted in the text.

## Data

All input data are public. Fetch scripts will be provided here; raw data is not committed.
Sources, as described in the manuscript:

- HamSCI amateur radio spots (RBN, PSKReporter, WSPRNet) from CEDAR Madrigal
- AE and SYM-H geomagnetic indices, and F10.7 solar radio flux, from NASA OMNI
- The Polar Vortex Index of Frissell et al. (2016), computed from MERRA-2 reanalysis
- Autoscaled ionosonde parameters from mid-latitude US stations of the Global Ionosphere Radio
  Observatory (GIRO)

## AI use

This repository is developed with AI assistance (Claude, Anthropic) under the HamSCI Generative
AI Use Agreement. Every AI-assisted session is logged in [`ai/ai_usage_log.md`](ai/ai_usage_log.md),
and AI-assisted commits are prefixed `[AI-assisted]`. The authors review and are responsible for
all code and results.

## Acknowledgments

This work is supported by the US National Science Foundation under Grant AGS-2045755. Further
support is acknowledged in the paper.

## License

MIT. See [`LICENSE`](LICENSE).
