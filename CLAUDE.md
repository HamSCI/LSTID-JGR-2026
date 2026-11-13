# LSTID-JGR-2026: Figure Code for the LSTID Climatology Paper

## Project Overview
This **public** repository holds the code that produces the figures and derived results of the
manuscript *Multi-Year Climatology of Large-Scale Traveling Ionospheric Disturbances Observed
with High Frequency Amateur Radio Receiving Networks Using a Novel Automated Detection
Algorithm* (Sanchez et al.), in preparation for *JGR: Space Physics*. It will be archived on
Zenodo when the paper is submitted, and the paper's Open Research section will cite that
archive. Its audience is reviewers and readers who want to reproduce the paper's figures, and
the authors themselves.

The LSTID detection algorithm is **not** in this repository. It is the separately released
[`HamSCI/hamsci_LSTID_detection`](https://github.com/HamSCI/hamsci_LSTID_detection) (v1.0,
Zenodo DOI in its README), which this repository uses as a pinned dependency. This repository
starts from that pipeline's output and the external index data, and ends at the paper's figures.

Planning, notes, and coordination live in the private orchestration repository
`HamSCI/LSTID-PRIVATE`, where this repository is the `code/` submodule. **Nothing private goes
here**: no unpublished collaborator data, no review correspondence, no draft manuscript text.

**Lead**: Diego F. Sanchez (KD2RLM), The University of Scranton (M.S. Software Engineering 2026)
**PI and advisor**: Nathaniel A. Frissell (W2NAF), The University of Scranton
**Collaborators**: the manuscript's co-authors, as listed in the manuscript
**HamSCI Working Group**: none
**Institution**: The University of Scranton
**Funder**: NSF AGS-2045755 (CAREER), primary. NASA support (80NSSC21K1772, 80NSSC23K0848) is
acknowledged in the detection code; whether it is acknowledged here is the authors' decision.
**Project period**: through submission and review of the paper (2026 to 2027)

## Project Goal
Anyone can regenerate every figure in the paper from public data and this code, with one
documented command per figure, and the version cited by the paper is archived on Zenodo.

## Standing Rules

@.claude/rules/ai-governance.md
@.claude/rules/hamsci-data.md

`.claude/rules/python-code.md` applies to `.py`, `pyproject.toml`, and `requirements*.txt`
through its own `paths:` frontmatter. This repository has no LaTeX; the LaTeX rule file is not
adopted.

## Figures

The paper's figures, in order of appearance in the manuscript (`main.tex` at `kd2rlm` commit
`3015b41`, 2026-11-13). The `\label` is the stable name; numbers shift as the draft changes.

| Label | File in the manuscript repo |
|---|---|
| `fig:ham_obs` | `plots/fig_ham_obs.jpg` |
| `fig:Pre_Processed` | `plots/fig_preprocessed_v3.png` |
| `fig:Sinfit` | `plots/fig_sinfit.png` |
| `fig:climatology_time_series` | `results_plots/Daily Sinusoidal Fit Parameters.png` |
| `fig:correlation_panel_amp_pv` | `results_plots/fig_correlation_panel_amp_pv.png` |
| `fig:seasonal_phase` | `results_plots/fig_seasonal_phase_amp_pv_ae_symh.png` |
| `fig:manual_auto_validation_2017` | `results_plots/fig_rmse_thresh.png` |
| `fig:monthly_lstid_vs_calm` | `results_plots/20150101-20221231_monthly_event_counts_combined.png` |
| `fig:event_probability_modulation` | `results_plots/20150101-20221231_event_probability_cpdfs_vertical.png` |
| `fig:correlation_panel_minskip_flux` | `results_plots/fig_correlation_panel_minskip_flux.png` |
| `fig:correlation_panel_minskip_fof2_es_summer_winter` | `results_plots/fig_correlation_panel_minskip_fof2_es_summer_winter.png` |

Keep this table in step with the manuscript. Each row should eventually name the script that
makes the figure and the inputs it reads.

## Working Conventions

**Reproducibility is the deliverable.** Every figure is made by a committed script from inputs
that are either public (fetched by a committed script) or produced by a committed step. No
figure depends on a notebook state, a local path, or a hand edit. Pin dependency versions,
including `hamsci_LSTID_detection` at its released tag.

**Data stays out of git.** Madrigal HDF5 files, pipeline caches, and index downloads are
fetched or regenerated. Small derived tables the figures read directly may be committed if the
authors decide so, with their provenance stated. Generated figures are not committed here; the
paper's copies live in the manuscript repository.

**Cite data and software.** HamSCI amateur radio spots via CEDAR Madrigal, OMNI (AE, SYM-H,
F10.7), MERRA-2 behind the Polar Vortex Index, GIRO ionosonde data, and the detection software,
each with access dates. See `.claude/rules/hamsci-data.md`.

**Commits.** Use `/commit`. Prefix AI-assisted commits with `[AI-assisted]`. Every change goes
on a feature branch and through a pull request; a human maintainer reviews and merges. Claude
never commits or pushes to `main` and never merges. Pushing a feature branch and opening its PR
is standing permission once the user approves the commit; any other push needs explicit
instruction. Never force-push or hard-reset.

**Releases are human decisions.** The Zenodo archive is made by a GitHub release with the
Zenodo integration enabled. The tag, the `CITATION.cff` author list, and the license are set by
the authors; Claude can prepare them but does not tag or release.

## Template Updates

This repository follows `HamSCI/ai_project_template`; `.claude/template.json` records the
template commit last adopted. A session-start hook checks for template changes once a day; when
it reports one, offer `/template-sync`.

## AI Governance

Every substantive AI session is logged in `ai/ai_usage_log.md` **before** the work is committed.
Use `/commit`, which enforces the ordering. Because this repository is public, its log is public
too, which is intended: it backs the AI-use disclosure in the paper and the Zenodo release.
