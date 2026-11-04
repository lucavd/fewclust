# cran-comments

## Submission

This is a new submission (fewclust 0.4.3).

## Test environments

- Windows Server 2022, R-devel (2026-11-03 r90835 ucrt), win-builder:
  0 errors, 0 warnings, 1 NOTE.
- Ubuntu 24.04 (noble), R 4.6.1: `R CMD check --as-cran`, with the CRAN
  incoming checks and their remote part switched on
  (`_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`): 0 errors, 0 warnings, 2 NOTEs.

## R CMD check results

- NOTE: New submission.
- NOTE (win-builder): possibly misspelled words in DESCRIPTION: DeRouen,
  Finos, Goeman, Graubard, Heagerty, Hemerik, Kauermann, Mancl, tibbles.
  These are the surnames of the authors of the cited methods and the name
  of the tibble data frames of the tidyverse; they are spelled correctly.
- NOTE (Ubuntu only): HTML validation of the manual was skipped because
  HTML Tidy is not installed on that machine; on win-builder the HTML
  manual check is OK.

## Notes for the reviewers

- The marginalized multilevel model of Heagerty (1999) is implemented in the
  package because the only other R implementation (MMLB) is not on CRAN.
  The test suite compares the fits with reference values computed by MMLB
  and stored as small fixtures.
- The small-sample sandwich corrections are implemented in the package.
  For the Kauermann-Carroll correction, the implementation in 'geesmv'
  applies a symmetric matrix square root to a matrix that is not symmetric;
  the package uses a similarity transformation that makes it symmetric first.
- The defaults follow a simulation study of these methods (Vedovelli et
  al., manuscript); the package's vignette summarises its findings.
