# AI Usage Log: LSTID-JGR-2026 (figure code)

This log records every substantive AI-assisted session on the project "LSTID-JGR-2026: figure code for the LSTID climatology paper (Sanchez et al.)".

Required by the HamSCI Generative AI Use Agreement, and by any institutional or funder policy
that applies to this project (see `.claude/rules/ai-governance.md`).

**This log is the source of truth for what AI did on this project.** A disclosure paragraph in a
manuscript, poster, or software release summarizes this log truthfully. Keep it current: an
entry written from memory weeks later is not a record.

## Entry format

```
## [YYYY-MM-DD HH:MM TZ]
- **Tool**: Claude (Anthropic), <exact-model-id>
- **Session Purpose**: What the session set out to accomplish
- **Sections/Files Affected**: Specific files, sections, or documents touched
- **Nature of Contribution**: Draft / Edit / Analysis / Code generation / Research / Scaffolding
- **Human Review Status**: Reviewed and verified / Partially reviewed / Pending review
- **Git Hash**: <filled in after committing>
```

The date and time come from the system clock via `date`, never from an estimate. The `/commit`
command produces this format and appends it in the right order.

---

<!-- Append new entries below this line, newest at the bottom. -->

## [2026-11-13 14:52 UTC]
- **Tool**: Claude (Anthropic), claude-opus-5-5
- **Session Purpose**: Turn the template copy into the public figure-code repository for Sanchez
  et al. (*JGR: Space Physics*), to be archived on Zenodo at submission; directed by
  N. A. Frissell (W2NAF).
- **Sections/Files Affected**: `CLAUDE.md` (project description, figure inventory from
  `main.tex` at kd2rlm `3015b41`, working rules), `README.md`, `CITATION.cff` (author list left
  for the authors), `.claude/rules/ai-governance.md` (Tier 2 University of Scranton; Tier 3 NSF,
  CAREER AGS-2045755, provisional NASA), `.claude/template.json` (project record, template
  `a09f960`), `ai/ai_usage_log.md` (template entries cleared), `.claude/rules/latex-writing.md`
  (removed; no LaTeX here)
- **Nature of Contribution**: Scaffolding, draft
- **Human Review Status**: Pending review
- **Git Hash**: d8155fb, branch `project-setup` (PR HamSCI/LSTID-JGR-2026#1)

## [2026-11-13 15:32 UTC]
- **Tool**: Claude (Anthropic), claude-opus-5-5
- **Session Purpose**: Correct the cited release of `hamsci_LSTID_detection`: the setup docs
  said "v1.0" (taken from a commit message), but the release is
  `V0.2.0_Sanchez_et_al_2026_MS` (Zenodo DOI 10.5281/zenodo.20186618), as the manuscript bib
  cites it. Verified with `gh release list`.
- **Sections/Files Affected**: `CLAUDE.md` (Project Overview), `README.md` (How this fits together)
- **Nature of Contribution**: Edit
- **Human Review Status**: Pending review
- **Git Hash**: 8b82594, branch `fix-detection-version` (PR HamSCI/LSTID-JGR-2026#2)

## [2026-11-17 05:50 UTC]
- **Tool**: Claude (Anthropic), claude-opus-5-5
- **Session Purpose**: Add the figure code for Figs 4–11 (moved from Diego's notebook) to the
  public repo so the figures can be regenerated: one script per figure plus `run_all.sh`, input
  location set by the `LSTID_DATA_DIR` environment variable (data kept out of git), a GIRO
  ionosonde download script, and docs (README, data/README, DATA_LICENSES, per-figure READMEs
  for Figs 1–3, which are made by harc_plot and hamsci_LSTID_detection).
- **Sections/Files Affected**: `2026_JGR_HamSCI_LSTID_Sanchez_figures/` (new): fig4–fig11
  scripts (Fig 4 output renamed to the manuscript's file name), `run_all.sh`,
  `requirements.txt` (adds scikit-learn for Fig 7), `.gitignore`, `README.md`,
  `DATA_LICENSES.md`, `data/README.md`, `data/gen_sza.py`, `data/download_giro.py`,
  `fig1/`–`fig3/README.md`
- **Nature of Contribution**: Code generation, edit (refactor of author's notebook code), documentation
- **Human Review Status**: Partially reviewed (all 8 scripts ran; Fig 10 r = −0.479 and −0.945
  match the paper; `download_giro.py` not yet tested against the GIRO server)
- **Git Hash**: 53be862, branch `add-figure-code` (PR HamSCI/LSTID-JGR-2026#3)

## [2026-11-23 04:43 UTC]
- **Tool**: Claude (Anthropic), claude-opus-5-5
- **Session Purpose**: Add analysis scripts that reproduce, from the real data, the derived
  numbers in the paper text that no figure script prints (tracked in the authors' review issues).
- **Sections/Files Affected**: `2026_JGR_HamSCI_LSTID_Sanchez_figures/analysis/` (new):
  `common.py`, `vortex_seasonal_tests.py`, `day_to_day_drivers.py`, `summer_sporadic_e.py`,
  `validation_2017_daily.py`, `qc_and_skip_windows.py`, `README.md`;
  `2026_JGR_HamSCI_LSTID_Sanchez_figures/run_all.sh` and `README.md` (run and document them)
- **Nature of Contribution**: Code generation, analysis, documentation
- **Human Review Status**: Partially reviewed (all five scripts run on the real data; numbers
  checked against the text; code not reviewed line by line)
- **Git Hash**: bec3f7e, 8a40e38, branch `add-analysis-scripts` (PR HamSCI/LSTID-JGR-2026#5)

## [2026-11-24 16:06 UTC]
- **Tool**: Claude (Anthropic), claude-opus-5-5
- **Session Purpose**: Remove three data limits that the figure scripts applied but the paper
  does not describe (amplitude 0-90 km, minimum skip distance 0-1600 km, F10.7 <= 200), so
  every figure uses only the quality control in paper Section 3.5 (R^2 >= 0.35, period
  1-4.5 hr). The F10.7 change keeps the removal of OMNI's 999.9 fill value.
- **Sections/Files Affected**: `fig4`-`fig11` scripts; `analysis/common.py`,
  `analysis/qc_and_skip_windows.py`, `analysis/vortex_seasonal_tests.py`
- **Nature of Contribution**: Code edit, analysis (all figure and analysis scripts rerun with and
  without the limits; outputs compared line by line)
- **Human Review Status**: Partially reviewed (removal requested by Diego after reviewing the
  before/after comparison; code diff not yet reviewed)
- **Git Hash**: [fill in after committing]
