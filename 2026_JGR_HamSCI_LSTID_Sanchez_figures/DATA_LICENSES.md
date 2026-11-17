# Data sources, licenses, and required acknowledgments

| Data | Source | Terms |
|---|---|---|
| `20140101-20230901_merged.csv` | This study, derived from the HamSCI Madrigal amateur radio dataset (CEDAR Madrigal, MIT Haystack Observatory) | Same license as this repository (TODO: author) |
| `2017_data_full_V2.csv` | This study (manual classification) | Same license as this repository (TODO: author) |
| `sym_ae.lst`, `omni2_daily_3hWFoJxMID.lst` | NASA GSFC Space Physics Data Facility, OMNIWeb (King & Papitashvili, 2005) | Public NASA data; acknowledge OMNIWeb/SPDF |
| `merra2_pvIndex.csv` | Derived from NASA GMAO MERRA-2 reanalysis by co-author V. L. Harvey | TODO: confirm with V. L. Harvey that the index may be redistributed; acknowledge GMAO MERRA-2 |
| `sza_daily.csv` | Computed by this study (`gen_sza.py`) | Same license as this repository |
| `giro/*.txt` | Lowell GIRO Data Center (DIDBase), Global Ionosphere Radio Observatory | CC BY-NC-SA 4.0 (see below) |

## GIRO ionosonde data

The GIRO files are redistributed unmodified under the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0), as stated in each
file header. Their use is subject to the Lowell GIRO Data Center Rules of the Road:
https://ulcar.uml.edu/DIDBase/RulesOfTheRoadForDIDBase.htm

Each file header states that its station's data provider must be acknowledged. Stations included:

| URSI code | Station |
|---|---|
| BC840 | Boulder, Colorado |
| AU930 | Austin, Texas |
| EG931 | Eglin, Florida |
| AL945 | Alpena, Michigan |
| WP937 | Wallops Island, Virginia |
| MHJ45 | Millstone Hill, Massachusetts |

TODO (author): add the acknowledgment text for the Lowell GIRO Data Center and each station's data
provider, matching the paper's Acknowledgments section.
