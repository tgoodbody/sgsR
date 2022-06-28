* `enhanced` -  `calculate_allocation()` improved documentation for output dataframe to make attributes more clear.

* `fixed` - `calculate_representation()` will now not plot bar chart twice.

* `fixed` - `existing` samples with other attributes will now not break sampling using `sample_ahels() / sample_clhs()` if values are `NA`. Variables are also added back to the sample output

# sgsR 1.0.0

* First CRAN release.

# sgsR 0.1.6

* Added comprehensive tests using `testthat` for most functions and `covr` reporting.

  - Resulting from tests: miscellaneous adjustments to many functions including small non-critical bug fixed, error and message improvements.
  
* Changed `forceSamp` to `force` in `sample_systematic()`.

* Removed `details` from `sample_coob()`.

* Improved ability to use `data.frames()` of samples as inputs for numerous algorithms.

* Improve consistency in error messages and `call. = FALSE` across the package.

* Added `existingna.shp` for example data where point are in `NA` locations.

# sgsR 0.1.5



* `strat_map()` can now map categorical srasters (gave an error before). Now also returns the categories associated with the categorical rasters in the lookup table with `details = TRUE`.

* Altered error handling for raster masking to be within `mask_access()` instead of individual sampling functions.

* Changes to `terra::distance()` & `terra::classify()` required slight modifications to `calculate_distance()` and `strat_breaks()`

* `sample_strat`
  - `strat_rule1()` and `strat_rule2()` functions to make code more concise. Fixed issue with `mindist` not setting distances between strata.

* `sample_ahels()`
  - Added `tolerance` parameter to `sample_ahels()` to allow users to define a tolerance around the desired sampling ratio (max `tolerance =  0.1`). This will allow the user to say "I ideally want the ratio to be `XX` but I'm OK if the ratio is `XX - threshold` if it means I don't need to add more samples".
  - Fix issue where an error occurred when `existing` was a data frame. 
  
* `sample_ahels` made `ahels_nSamp()` and `ahels_threshold()` functions to make code more concise.

* Altered how `existing` (crosses) and `new` (circles) samples are plotted in `sample_strat()` and `sample_ahels()`.

* Added internal utility functions `allocate_prop / allocate_optim / allocate_manual / allocate_equal / allocate_existing / allocate_force` and integrated them into `calculate_allocation()` to make code more succinct and purposeful.

# sgsR 0.1.4

* Fixed issue in `calculate_allocation()` where too many samples would be allocated (compared to the user-defined `nSamp`) due to using `ceiling()` instead of `round()` during proportional and optimal allocation.

* Added `allocation = "manual"` to `calculate_allocation()`. The parameter `weights` was added (mandatory for `allocation = "manual"`), where users can provide a numeric vector of relative weightings to strata. `sum(weights)` must equal 1.

* Added `weights` parameter to `sample_strat()` to allow for `"manual"` allocation.

* Allow `buff_inner` to be `NULL` when providing `access` to isolate samples. This allows users to define only a maximum distance (should they wish to) that samples can be from `access` but not specify a minimum distance.

# sgsR 0.1.3

* Updates names for internal package data in `inst/exdata` and corresponding examples and vignettes.

* Added contingencies for `sample_ahels`, `sample_strat`, `calculate_allocation`, `calculate_coobs` to
allow existing samples to fall in areas of `NA` and not cause algorithms to fail or bug.

* `extract_metrics` & `extract_strata` produce more robust error messages & will now generate a message if 
existing samples are co-located with strata or metric values that are `NA`.

* `sample_clhs` updated to fix issue where existing samples would not be appended to output.

* `strat_quantiles` now provides a stratum look-up table when `details = TRUE` to allow user to see exact metric break
points.

* Updated documentation for many algorithms to fix grammatical and spelling errors.

# sgsR 0.1.2

* Removed dependencies for `stringr`, `tidyselect`, `rlang`, `RColorBrewer`, `magrittr`, `glue`, `Rcpp`.

* Improved documentation examples to have read speeds in line with <5 seconds for CRAN submission.

* Improve processing speed for `sample_coobs()` and `utils-mat.R` courtesy of Jean-Romain Roussel.

* Added example subset file to `inst/extdata` to improve example processing speeds for certain algorithms.

# sgsR 0.1.1

* Added a `NEWS.md` file to track changes to the package.
