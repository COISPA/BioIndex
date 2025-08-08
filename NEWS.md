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
