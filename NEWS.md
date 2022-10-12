# sgsR 1.3.1

`fixed` - `sample_systematic()` added random translation to sampling grid to ensure unbiased sampling.

# sgsR 1.3.0

`new` - `sample_existing()` has been added. This algorithm sub-samples and `existing` sample using internal latin hypercube sampling. Constraints in the form of the `cost` parameter akin to `sample_clhs()` exist. Sub-sampling can be performed on `existing` samples alone, or using population level `mraster` distributions.

`new` - `mask_existing()` - internal function for masking `exising` samples using `access` buffers.

# sgsR 1.2.1

`fixed` - `sample_systematic()` bug where `cellsize` values that caused no samples to intersect with the raster would cause `extract_metrics()` to provide an error about existing not having any samples. Added a check for intersection and corresponding error message.

`enhancement` - Added `quiet` to `extract_metrics()` & `extract_strata()` to allow internal use without messages.

# sgsR 1.2.0

* `fixed` - `strat_kmeans()` bug related to `terra` where the re-assignment of values to the output raster was causing issues. R Hijmans kindly suggested the edit made.

* `fixed` - `sample_ahels()` - bug where extra attributes in `existing` would cause the algorithm to crash when re-merging after sampling.

* `fixed` - `strat_quantiles()` no longer plots histogram / scatter plot when using `plot = TRUE`. Now correctly adds this to details list when `details = TRUE`.

* `new sampling method` - Added `sample_nc()` based on the algorithm described in [Melville & Stone (2016)](https://doi.org/10.1080/00049158.2016.1218265)

    * This algorithm uses kmeans clustering where the number of clusters is equal to the desired number of samples. Cluster centers are located, which then prompts the nearest neighbour raster pixel for each cluster to be located (assuming default `k` parameter). These nearest neighbours are the output samples. Visualization of the centers and samples can be dispayed if `details = TRUE` is used and `$kplot` is plotted.

* `fixed` - `sample_systematic()` now has inherent randomness for lower left corner of the tessellation.

* `fixed` - `strat_kmeans()` solved issue where only first raster layer was being involved in stratification.

# sgsR 1.1.0

* `enhanced` - `sample_strat()` - added parameter `method` that allows users to choose between `"Queinnec"` (default method implemented in previous sgsR versions) and `"random"` (stratified random sampling). The random method ignores much of the functionality of the algorithm to allow users to use standard stratified random sampling approaches without the use of a focal window to locate contiguous stratum cells.

* `fixed` - `sample_strat()` factor handling improvement - GitHub issue #18

* `enhanced` -  `calculate_allocation()` improved documentation for output data frame to make attributes more clear.

* `fixed` - `calculate_representation()` will now not plot bar chart twice & `NA` values in existing will not be removed.

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
