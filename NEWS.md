# BioIndex v0.4.04
  * added IUT functions
# BioIndex v0.5.01
  * added functions to perform analyses on TE table (ALK and LW)
# BioIndex v0.5.02
  * improved merge.TATBTC function
  * improved outputs saving in MIW function
  * improved outputs saving in LFD, Lquant functions
  * improved outputs saving in index_on_grid, sex_ratio_on_grid functions
  
# BioIndex v0.6.00
  __Fixes for RBDFIS III__
  
  * Optimised the fiction `MIW()` for saving and outputing the plots in the correct way in the console.
  * function `MIW()`corrected the functioning of verbose and save parameters
  * Improved `bubbleplot_RS_by_hauls()` to allow offline plotting: if `res = NA`, the function now uses an internal dataset `med_bathy` (class `"bathy"`) with preloaded bathymetry for the Mediterranean and Black Sea (0 to -1000 m). Online data retrieval using `marmap::getNOAA.bathy()` is still available by specifying a resolution value (e.g., `res = 1`). This dual approach ensures reproducibility and fast plotting in both offline and online environments, reducing the computational times accessing to the embedded dataset.
  * Improved `bubble_plot_by_haul_indexes()` to allow offline plotting: if `res = NA`, the function now uses an internal dataset `med_bathy` (class `"bathy"`) with preloaded bathymetry for the Mediterranean and Black Sea (0 to -1000 m). Online data retrieval using `marmap::getNOAA.bathy()` is still available by specifying a resolution value (e.g., `res = 1`). This dual approach ensures reproducibility and fast plotting in both offline and online environments, reducing the computational times accessing to the embedded dataset.
  * Improved `hauls_position()` to allow offline plotting: if `res = NA`, the function now uses an internal dataset `med_bathy` (class `"bathy"`) with preloaded bathymetry for the Mediterranean and Black Sea (0 to -1000 m). Online data retrieval using `marmap::getNOAA.bathy()` is still available by specifying a resolution value (e.g., `res = 1`). This dual approach ensures reproducibility and fast plotting in both offline and online environments, reducing the computational times accessing to the embedded dataset.
  * function `bubbleplot_RS_by_hauls()` corrected the functioning of verbose and save parameters.
  * function `bubble_plot_by_haul_indexes()` corrected the functioning of verbose and save parameters.
  * function `hauls_position()` corrected the functioning of verbose and save parameters.
  * function `LFD()` added verbose messages to enhance traceability and user feedback. Ensured consistent and explicit handling of missing data scenarios via informative verbose messages instead of cat() statements.
  * function `index_on_grid()` added verbose messages to enhance traceability and user feedback. Ensured consistent and explicit handling of missing data scenarios via informative verbose messages instead of cat() statements.
  * function `sex_ratio_on_grid()` added verbose messages to enhance traceability and user feedback. Ensured consistent and explicit handling of missing data scenarios via informative verbose messages instead of cat() statements.
  * function `overlayGrid()` added verbose messages to enhance traceability and user feedback. Ensured consistent and explicit handling of missing data scenarios via informative verbose messages instead of cat() statements.
  * Added standalone function `merge_TATB()` to merge TA and TB tables only,  
    including full RoME QC. Mirrors legacy behaviour but faster on large datasets thanks to vectorization of for loops.
  * Added standalone function `merge_TATC()` to merge TA and TC tables only,  
    including full RoME QC. Mirrors legacy behaviour but faster on large datasets thanks to vectorization of for loops.
  * Optimized `merge_TATBTC()` vectorising for loops to reduce run-time while keeping identical
    numerical results.
  * `BioIndex()` function modified to create a report filename that now includes the MEDITS code of the analyzed species, the refernece GSA and the lower and upper depth values (e.g., _BioIndex_results_ARISFOL_GSA18_Depth200-800m_2025-08-08_h09m15s24.zip_). This makes it easier to organize and identify files, especially when running analyses for multiple species or depths (according to outcomes of RDBFIS II second training).
  * `BioIndex()` documentation modified to specify the *mm units* required for both recruits and spawners cutoff thresholds (according to outcomes of RDBFIS II second training).
  * `bubbleplot_RS_by_hauls()` documentation modified to specify the mm units required for both recruits and spawners cutoff thresholds (according to outcomes of RDBFIS II second training).
  * `index_recr()` documentation modified to specify the mm units required for recruits cutoff threshold (according to outcomes of RDBFIS II second training).
  * `index_spawn()` documentation modified to specify the mm units required for spawners cutoff threshold (according to outcomes of RDBFIS II second training).

# BioIndex v0.6.01
  __Fixes for RBDFIS III__
  * `index_on_grid()` modified to auto select the bathymetry contour polygon to be used in the plots.

# BioIndex v0.6.02
  __Fixes for RBDFIS III__
  * The stratification tables relative to the rapa whelk beam trawl survey, conducted in the the Black Sea by Romania and Bulgaria, are included in the library: strata_scheme_rapana and stratification_rapana.

# BioIndex v0.6.03
Added an embedded Shiny application to the package: run_BioIndex_app(). It provides a graphical user interface for running BioIndex analyses directly from the package environment.

# BioIndex v0.6.4
  __CRAN Submission Release__
  
  * **Robustness and Bug Fixes**:
    * Fixed stratum variance and abundance/biomass calculations in `indices_ts.R`, `index_recr.r`, `index_spawn.r`, `index_ts_F.r`, `index_ts_M.r`, and `MIW.r` to robustly filter `NA` values, handle `NA`/`NaN` in `cv_strata` calculations, and prevent division by zero in GSAs with a single haul.
    * Corrected the sex ratio variance calculation in `sex_ratio.r` by removing the incorrect square root from the variance formula.
    * Added automated directory creation for the output folder in `Lquant.r` before saving files.
    * Corrected `@return` documentation in `bubbleplot_RS_by_hauls()` to describe the returned list of plots, initialized them to `NULL`, and returned them correctly.
    * Refactored `convert_coordinates.r` to safely skip `NA` coordinate records and avoid subscripted assignments error.
  
  * **CRAN Policy Compliance**:
    * Standardized directory handling: functions now default to `tempdir()` when `wd = NA`, ensuring compliance with CRAN's file system policies regarding unauthorized write access. Removed all `setwd()` calls from internal logic to use absolute paths (`file.path()`).
    * Refactored console output: replaced all `cat()` and `print()` calls with `message()` wrapped in `if (verbose)` checks to ensure a silent default behavior.
    * Improved state management: added `on.exit(par(oldpar))` and secured working directory restoration where applicable.
    * Documentation updates: added missing `@return` tags (e.g., in `hauls_position()`) as requested by CRAN reviewers.
    * Example Activation: unwrapped examples from `\donttest{}` in `aggregate_gsas()`, `merge_TATB()`, `merge_TATC()`, and `merge_TATBTC()` to allow automatic testing.
  * **Dependency Refinement**: Internalized data validation logic based on RoME (v0.2.3) to eliminate external dependencies on non-CRAN repositories.
  * **Bug Fixes and Stability**:
    * Fixed non-ASCII characters in `continent` dataset.
    * Refactored `class()` checks to use `inherits()` for robust type checking.
    * Added `globalVariables` declarations to resolve R CMD check notes.
    * Standardized function signatures to consistently include `wd` and `save` parameters.
    * Fixed function signature parameter mismatches:
      * Added missing `GSA` (default 10) and `wd` (default NA) arguments to the `IUT()` function signature.
      * Added missing `verbose` (default FALSE) argument to the `spearman()` function signature.
    * **Enhanced Data Consistency**:
        * Implemented explicit GSA and Country filtering in `BioIndex()` to ensure all downstream analysis functions use the correct subset of data and stratification metadata, preventing errors when multiple GSAs are present in input files.
        * Added robust validation for depth range boundaries in `indices_ts()`, `sex_ratio()`, `MIW()`, `LFD()`, and other core functions. The package now provides clear error messages when the selected depth range does not match the strata defined for the GSA.
        * Updated the embedded Shiny application to dynamically set depth range defaults (Depth lower/upper) based on the selected GSA and Country, preventing common user configuration errors.
        * Improved Shiny UI: replaced numeric inputs for "Depth lower" and "Depth upper" with dynamic dropdown menus (select inputs) populated directly from the stratification table boundaries.
        * Refined Shiny UI default logic: "Depth lines" default values are now intelligently extracted from the actual stratification cutoffs (selecting the two extremes and the most central cutoff) rather than using an arbitrary mathematical average, resulting in more meaningful default map contours.
        * **Data Integration and Reproducibility**:
            * Fixed a critical issue in `merge_TATC` and `merge_TATBTC` where decimal coordinates were not correctly assigned, causing failures in spatial aggregation routines (e.g., in `bubbleplot_RS_by_hauls`).
            * Enhanced stability of merge functions to robustly handle empty data subsets, preventing `cbind` mismatches.
            * Optimized example execution times: wrapped the heavy global workflow example `BioIndex()` in `\donttest{}` and subset individual mapping examples (`bubble_plot_by_haul_indexes()` and `hauls_position()`) to a single year (2016) to ensure swift execution times under CRAN checks.
            * Standardized all documentation examples to utilize internal GSA 10 data, ensuring reliable execution and full compliance with CRAN's automated checking environment.
        * **Graphics and UI Stability**:
            * Eliminated all explicit `print()` calls for `ggplot2` objects within analytical functions, relying exclusively on `ggsave()` for output. This prevents fatal 'figure margins too large' and 'pin' errors when functions are executed via the Shiny app on small active graphics devices.
            * Silenced rogue base R graphics rendering in functions like `Lquant()` that previously attempted to plot to the active screen device by default.
            * Refined graphical parameter restoration logic: `on.exit(suppressWarnings(par(oldpar)))` now explicitly strips the `new` attribute from `oldpar`, preventing 'calling par(new=TRUE) with no plot' warnings after execution.
            * Suppressed non-critical `ggplot2::stat_contour` warnings (such as 'Zero contours were generated') across spatial plotting functions to maintain a clean console output as per CRAN policies.
