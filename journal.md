# Journal Requirements

## Active journal

- Journal: Journal of Open Source Software (JOSS)
- Article type: Software paper
- Submission guide: https://joss.readthedocs.io/en/latest/submitting.html
- Paper format: https://joss.readthedocs.io/en/latest/paper.html
- Review criteria: https://joss.readthedocs.io/en/latest/review_criteria.html
- Review checklist: https://joss.readthedocs.io/en/latest/review_checklist.html
- Guide source: official JOSS documentation
- Last verified: `2026-10-30 09:44` (Asia/Shanghai)
- Verified by: Codex, from the official JOSS documentation
- Paper source: `paper.md`
- Bibliography: `paper.bib`
- Submission status: preparing the public review issue after release checks

## Hard requirements

### Scope and significance

- The software must be open source under an OSI-approved license.
- The repository must be publicly clonable and browsable without registration.
- The issue tracker must be public and accept new issues.
- The paper must describe research software, not new research results enabled by
  the software.
- The software must be feature-complete, packaged to R community standards, and
  maintained for reuse rather than one-off analysis.
- A meaningful contribution to the research community is required; a minor
  utility, thin API client, or single-function package is not eligible.
- Evidence of research use is required, at minimum by the developers. Private
  workflow use is acceptable when its evidence can be shown to the editor.

### Pre-review gates

- Public development history must exceed six months.
- Development must show ongoing iteration over time rather than a concentrated
  code dump.
- Single-author projects must show multiple open-source signals: sustained
  commits, tags or changelog, tests and CI, documentation, contribution
  guidance, and support or governance expectations.
- External issues, pull requests, and community engagement are positive signals
  but are not mandatory for a single-author project with other evidence.

### Paper format

- Length: 750-1750 words.
- Required metadata: title, tags, authors, affiliations, date, and bibliography.
- Required sections: Summary, Statement of need, State of the field, Software
  design, Research impact statement, AI usage disclosure, Acknowledgements, and
  References.
- The Summary must explain the software to a non-specialist audience.
- The State of the field must compare related tools and justify building rather
  than contributing to an existing project.
- The paper must not contain API documentation or a full software manual.
- References must use full venue names and include related software.
- Financial support and conflicts of interest must be acknowledged.

### AI usage disclosure

The paper must state the tools and models used, the scope of assistance, and how
the authors validated the outputs. Human authors remain responsible for all
technical and scholarly decisions.

### Review and archiving

- The paper is reviewed publicly in a GitHub issue by at least two reviewers.
- The authors must respond to reviewer comments within the requested review
  window.
- After successful review, create a tagged release and archive it with Zenodo,
  Figshare, or another accepted repository.
- Record the release version, archive DOI, and code availability in the review
  issue.

## Local compliance assets

- `paper.md`: JOSS manuscript with all required sections.
- `paper.bib`: verified bibliography.
- `CONTRIBUTING.md`: contribution and testing workflow.
- `CODE_OF_CONDUCT.md`: community behavior expectations.
- `SUPPORT.md`: support, issue, privacy, and maintenance guidance.
- `.github/ISSUE_TEMPLATE/`: structured bug and feature reports.
- `.github/workflows/draft-pdf.yml`: JOSS draft-PDF build.
- `.github/workflows/R-CMD-check.yaml`: multi-platform package checks.
- `JOSS_SUBMISSION.md`: pre-submission status and known editorial risk.

## Change log

- `2026-10-30 09:44` - Activated JOSS and archived the previous Bioinformatics
  Application Note requirements.
- `2026-10-28 21:47` - Initial journal record from the official OUP author
  guidelines.

## Archived journal requirements

### Bioinformatics Application Note (archived 2026-10-30)

- Journal: Bioinformatics (Oxford University Press)
- Article type: Application Note
- Guide URL: https://academic.oup.com/bioinformatics/pages/author-guidelines
- Guide verified: `2026-10-28 21:47` (Asia/Shanghai)
- Main text: up to four pages, approximately 2,600 words, or 2,000 words plus
  one figure.
- Abstract: structured as Summary, Availability and Implementation, Contact,
  and Supplementary information; recommended maximum 150 words.
- Scope: trivial utilities are not eligible, and installation must not require
  a significant investment of time.
- Figures: vector formats preferred; raster at 350 dpi or higher and line art at
  1200 dpi.
- Required declarations: data availability, code availability, funding,
  competing interests, and CRediT.
- Availability: software must remain available for at least two years.
- Fallback target at that time: Bioinformatics Advances.
