## Resubmission

This is a resubmission. All modifications requested by CRAN reviewers across the different review rounds have been fully implemented. The package has undergone multiple revision cycles addressing policy compliance, documentation completeness, dependency management, and code robustness; every individual point raised has been resolved, as detailed below.

*   **Dependency on non-CRAN repository removed**: The 'RoME' package has been removed from the `Suggests` field. Data validation logic previously delegated to 'RoME' has been fully internalized.

*   **`\donttest{}` usage justified**: The four spatial/cartographic functions (`BioIndex()`, `bubble_plot_by_haul_indexes()`, `hauls_position()`, `bubbleplot_RS_by_hauls()`) remain wrapped in `\donttest{}` because their execution time exceeds CRAN thresholds (> 10s on Windows, > 5s on Debian) due to `ggplot2` bathymetric contour rendering. All other functions have active, executable examples.

*   **Commented-out example code resolved**: The commented example in `overlayGrid.Rd` has been uncommented and is now fully executed during checks.

*   **`@return` documentation completed**: Missing `@return` tags have been added to all functions that were flagged, including `hauls_position()` and `bubbleplot_RS_by_hauls()`.

*   **Non-ASCII characters fixed**: Non-ASCII characters in the `continent` dataset have been resolved.

*   **`class()` checks refactored**: All `class()` equality checks have been replaced with `inherits()` for robust type checking.

*   **Namespace conflict resolved**: The conflict between `terra::extract()` and `magrittr::extract()` has been resolved.

*   **CRAN file system policy**: All functions now default to `tempdir()` when `wd = NA`. All `setwd()` calls have been removed from internal logic.

*   **Console output silenced**: All `cat()` and `print()` calls have been replaced with `message()` wrapped in `if (verbose)` guards, ensuring silent default behavior.

*   **Graphics state restored**: `on.exit(par(oldpar))` has been added wherever graphical parameters are modified.

## Test environments

- Local Windows 10, R 4.5.1
- win-builder (R-devel, 2026-05-23 r90071): `devtools::check_win_devel()`

## R CMD check results

### Local check

0 errors | 0 warnings | 0 notes

### win-builder (R-devel)

0 errors | 0 warnings | 1 note

The single NOTE reported by win-builder is:

```
checking CRAN incoming feasibility ... NOTE
New submission

Possibly misspelled words in DESCRIPTION:
  GSAs (26:29)
  MEDITS (3:46, 16:40)
  RDBFIS (23:20)
```

These are not misspellings. They are established technical acronyms in the domain of fisheries science:
- **MEDITS**: MEDIterranean Trawl Survey — the name of the international trawl survey program this package supports.
- **GSAs**: Geographical Sub-Areas — the standard FAO spatial units used in Mediterranean fisheries management.
- **RDBFIS**: Regional DataBase and FIsheries Information System — the EU data collection framework.

All three terms have been added to the package `WORDLIST` file.
