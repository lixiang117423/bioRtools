# JOSS submission status

Official guidance checked on 30 October 2026:

- Submission guide: <https://joss.readthedocs.io/en/latest/submitting.html>
- Paper format: <https://joss.readthedocs.io/en/latest/paper.html>
- Review criteria: <https://joss.readthedocs.io/en/latest/review_criteria.html>
- Review checklist: <https://joss.readthedocs.io/en/latest/review_checklist.html>

## Required files

- [x] `paper.md` in the repository root with JOSS YAML metadata.
- [x] `paper.bib` with verified bibliography entries.
- [x] OSI-approved MIT text in `LICENSE.md` and `COPYING`.
- [x] Installation, usage, API, and workflow documentation.
- [x] Automated tests and continuous integration.
- [x] Contribution, support, and community-conduct guidance.
- [x] Draft-paper GitHub Actions workflow.

## Pre-review evidence

- Public repository created: 28 July 2025.
- Development history: more than 300 commits over 15 months.
- Version history: more than 15 tags and a maintained changelog.
- Automated tests: testthat 3e suite with public-data end-to-end workflows.
- Developer adoption: 944 R, R Markdown, and Quarto files across 22 active
  research projects call exported `bioRtools` functions.
- Documentation: README, roxygen reference, pkgdown site, and workflow vignettes.
- Test gate: 0 errors, 0 failures, 0 warnings; 6 tests skip because the optional
  `gggenes` package is unavailable.
- Package check: `R CMD check --no-manual` completed with `Status: OK`; missing
  optional Suggests produced INFO entries only.
- JOSS paper: the official `openjournals/inara` container compiled a 4-page
  draft from `paper.md` and `paper.bib` without a build failure.
- pkgdown: `pkgdown::check_pkgdown()` reported no problems.

## Before opening the JOSS review issue

1. Push the reviewed default branch and confirm GitHub Actions pass.
2. Create the `v2.0.0` GitHub release from the tested commit.
3. Confirm the pkgdown Pages deployment is active.
4. Provide the editor with a concise description of private-workflow research
   use, because the paper reports developer adoption without claiming external
   citations.
5. Archive the accepted release in Zenodo or Figshare and add the DOI after
   JOSS review acceptance; JOSS requests this archive before final acceptance.

## Known editorial risk

The project is currently single-author and has no external peer-reviewed
citations. JOSS permits single-author submissions and private-workflow adoption,
but the editor may require stronger public evidence of research use. The paper
states this limitation explicitly rather than inventing adoption.
