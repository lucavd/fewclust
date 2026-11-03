# cran-comments

## Submission

This is a new submission (fewclust 0.4.0).

## Test environments

- Ubuntu 24.04 (noble), R 4.6.1: `R CMD check --as-cran` -- 0 errors, 0 warnings, 1 NOTE.
  The NOTE says that the HTML validation of the manual was skipped because
  HTML Tidy is not installed on the test machine.  The remote checks of the
  incoming feasibility step (URLs) were switched off, as the repository is
  made public with the submission.

## Notes for the reviewers

- The `URL` and `BugReports` fields point at the public repository that
  hosts the package; it is made public together with this submission, so
  the incoming check may report the two URLs as not yet reachable.
- The marginalized multilevel model is implemented inside the package
  because the only other implementation (MMLB) is not on CRAN.  The test
  suite compares the fit with reference values computed by MMLB and
  stored as small fixtures.
- The small-sample sandwich corrections are implemented inside the package
  because the pre-existing implementation on CRAN was found numerically
  unstable at small samples in the simulation study that motivates this
  package, and computes a different quantity for the Kauermann-Carroll
  correction.
