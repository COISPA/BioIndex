## Test environments

- Local Windows 10, R 4.5.1
- win-builder (R-devel): `devtools::check_win_devel()`

## R CMD check results

0 errors | 0 warnings | 1 note

- NOTE: checking top-level files ... NOTE
  Non-standard file/directory found at top level: 'cran-comments.md'
  (This file is now added to .Rbuildignore and should not appear in the final tarball).

## Submission summary

This is a new submission for the BioIndex package.
BioIndex is designed to support the standardized analysis of MEDITS trawl survey data and the calculation of biological indicators.

### Key Changes in v0.6.04:
- Resolved `DESCRIPTION` field policy: Description now starts with a neutral verb instead of the package name.
- Fixed non-ASCII characters in `continent` dataset.
- Refactored `class()` checks to `inherits()` for robust type checking.
- Resolved namespace conflict for the `extract()` function between `terra` and `magrittr`.
- Updated author metadata and standardized documentation.
- Added `WORDLIST` and `inst/WORDLIST` to address spell-check warnings on technical acronyms (MEDITS, GSAs, RDBFIS, etc.).
- Fixed deprecated `ggplot2` arguments (`size` replaced by `linewidth`).
- Wrapped the main function example in `\dontrun{}` to ensure it is never executed during automated checks, avoiding failures related to internet connectivity (NOAA downloads) and execution time limits.

### Dependencies:
- **RoME**: BioIndex depends on the RoME package, which is currently undergoing a separate new submission process.

## Misspelled words (Technical Acronyms)
The following words reported as potentially misspelled are technical acronyms or domain-specific terms and have been added to the WORDLIST:
- BioIndex
- GSAs
- MEDITS
- RDBFIS
