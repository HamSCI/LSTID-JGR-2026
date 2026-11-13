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
- **Git Hash**: [fill in after committing]
