# sgsR 1.4.4

`fixed` - `sample_strat()` - Was not taking samples from other strata into account when calculating `mindist` between sample units. This has now been corrected for both `Queinnec` and `random` methods. From Tommasso Trotto #33

`enhanced` - `calculate_distance()` - added `slope` parameter to allow for calculation of slope distance. From Nic #32

# sgsR 1.4.3

`enhanced` - `sample_ahels()` - If attributes in `existing` contain numeric data -- and match columns in `mraster` -- these data will be used instead of extracting from `mraster`.

# sgsR 1.4.2

`added` - `utils-write` - made writing functions for samples, rasters, and tables to make code more succinct.

`enhanced` - `sample_existing` - A major change to `sample_existing()` has been implemented. Prior to `sgsR` v1.4.2, `sample_existing()` only leveraged `sample_clhs()` functionality. New functionality has been added to allow users to define the sub-sampling method of their choice out of (`sample_srs()`, `sample_balanced()`, `sample_strat()`, `sample_clhs()`). To do so, additional internal functions (`sample_existing_strat()`, `sample_existing_srs()`, `sample_existing_balanced()`, `sample_existing_srs()`) have been added that take in `existing` data and perform sampling based on those data. For `sample_existing_strat()`, new allocation algorithms that take `existing` as inputs needed to be written that are found in `utils-allocation-existing()`. Additionally, more utility functions that check `existing` data, and prepare `existing` data for sub-sampling have been developed and implemented in `utils-existing()`. Unit tests for all functions have been added. The `plot` parameter has been removed from `sample_existing()` for now. Likely to be added again later.

`enhanced` - Sampling Vignette - Added content about `sample_existing()`

`fixed` - `strat_map()` - Fixed improper use of `terra::hasValues`. From Robert Hijmans #31

# sgsR 1.4.1

`fixed` - `sample_ahels()` - gave an erroneous error when `matrices` was provided and `nQuant` did not match. Changed to allow only `matrices` to be provided.

`fixed` - `extract_strata()` - added error message if `sraster` had multiple layers #28

`fixed` - `extract_metrics()` - removed code where `mraster` layers with `ID` as name were removed. From Tommasso Trotto #27

`enhanced` - `calculate_representation()` - Added parameter `drop = NULL`. If `!=NULL` then filtering is applied where strata with under drop frequency are dropped from resulting dataframe and plot.  #28

# sgsR 1.4.0

`added` - `sample_sys_strat()` - Systematic stratified sampling. Using same functionality as `sample_systematic()` but takes an `sraster` as input and performs sampling on each stratum iteratively.

`enhanced` - `strat_breaks()` - Vectorized function to allow for any number of input mraster layers and a corresponding number of breaks vectors (as list in respective order as `mraster` layers). Removed `mraster2` & `breaks2`. Users can now supply an `mraster` with as many layers as they wish. Added `map` to allow for creating a combined (mapped) stratification output (`strata`). Internal function `calculate_breaks()` was added that facilitates vectorization.

`enhanced` - `strat_quantiles()` - Vectorized function to allow for any number of input mraster layers and a corresponding number of breaks vectors (as list in respective order as `mraster` layers). Removed `mraster2` & `nStrata2`. Users can now supply an `mraster` with as many layers as they wish along side `nStrata` as a list with `length(nStrata) == terra::nlyr(mraster)`. `nStrata` can be either a scalar integer representing the number of desired output strata, or a numeric vector of probabilities between 0-1 demarcating quantile break points. The `nStrata` list can be a mix of these (e.g. `nStrata = list(c(0.1,0.8,1), 4, 9)` where `mraster` would have 3 layers) to allow users to define both explicit quantile breaks or a desired strata number that is converted to quantiles breaks internally. Added `map` to allow for creating a combined (mapped) stratification output (`strata`). Internal functions `calculate_quantile_breaks() / quantile_breaks_integer() / quantile_breaks()` were added that facilitate vectorization.

`enhanced` - `strat_map()` - Vectorized function to allow for any number of input mraster layers and a corresponding number of breaks vectors (as list in respective order as `mraster` layers). Removed `raster2`. Users can now supply an `sraster` with as many layers as they wish. Thank you, Tommaso Trotto.

`enhanced` - Updated vignettes and documentation to account for vectorized functionality of the above functions.

`enhanced` - `plot_scatter()` - now visualizes with `viridis` colour scheme.

`fixed` - Added new citation information for upcoming manuscript release.

`fixed` - `extract_strata() / extract_metrics()` - CRS for `existing` will now be maintained when provided as an `sf` and the `mraster` CRS will be assigned otherwise. In addition, all sampling functions will maintain CRS from `existing` if possible, otherwise CRS from `sraster/mraster` are used for output samples.

# sgsR 1.3.4

`fixed` - `sample_srs() / sample_strat(method = "random")` - First sample unit was always duplicated From Tommaso Trotto.

`added` - `plot_scatter()` - Internal function. Scatter plot visualizing relationship between 2 `mraster` metrics with `existing` samples superimposed. 

# sgsR 1.3.33

`fixed` - `sample_strat()` - `srasters` with categorical values were crashing the algorithm due to inability to combine facter and non-factor values from Evan Muise.

`fixed` - `sample_systematic()` - Added more contingency for `cellsize` values that resulted in empty sample output.

# sgsR 1.3.32

`fixed` - `strat_map()` - `stratamapped` was outputting as character and not as integer or character depending on input strata type as intended from Tommaso Trotto.

# sgsR 1.3.31

`enchancement` - `calculate_pcomp()` - Added `maxcells` parameter based on suggestion from R. Hijmans.

`fixed` - `sample_systematic()` - Fixed issue related to ATLAS Blas and CRAN errors with suggestions and support from R. Hijmans.

# sgsR 1.3.3

`fixed` - CRAN issue where errors were encountered when run on ATLAS instances.

`enhancement` - Edited vignettes and documentation for clarity.

`enhancement` - `sample_srs()` - Added message to tell users when `nSamp` sample units were unable to be allocated. From Evan Muise.

# sgsR 1.3.21

`fixed` - `strat_quantiles() / strat_kmeans()` - solved issue where correct number of strata & float strata values were being output.

`enhancement` - `sample_existing()` - made it so extra attributes are passed to output when `raster` is provided. Added additional unit tests and updated documentation.

# sgsR 1.3.2

`enhancement` - `README.Rmd` and vignettes have been updated.

`enhancement` - `sample_systematic()` - changed how tessellation was used internally and visualized during plotting.

`fixed` - `strat_map()` - #20 greatly simplified algorithm using suggestion from R. Hijmans (added as author to algorithm). Issue was related to level matching with categorical variables.

`fixed` - `sample_ahels()` - swapped plot marker icons to be consistent with other functions.

# sgsR 1.3.1

`fixed` - `sample_systematic()` - added random translation to sampling grid to ensure unbiased sampling.

# sgsR 1.3.0

`new` - `sample_existing()` - has been added. This algorithm sub-samples and `existing` sample using internal latin hypercube sampling. Constraints in the form of the `cost` parameter akin to `sample_clhs()` exist. Sub-sampling can be performed on `existing` samples alone, or using population level `mraster` distributions.

`new` - `mask_existing()` - internal function for masking `exising` samples using `access` buffers.

# sgsR 1.2.1

`fixed` - `sample_systematic()` - bug where `cellsize` values that caused no samples to intersect with the raster would cause `extract_metrics()` to provide an error about existing not having any samples. Added a check for intersection and corresponding error message.

`enhancement` - Added `quiet` to `extract_metrics()` & `extract_strata()` to allow internal use without messages.

# sgsR 1.2.0

* `fixed` - `strat_kmeans()` bug related to `terra` where the re-assignment of values to the output raster was causing issues. R Hijmans kindly suggested the edit made.

* `fixed` - `sample_ahels()` - bug where extra attributes in `existing` would cause the algorithm to crash when re-merging after sampling.

* `fixed` - `strat_quantiles()` - no longer plots histogram / scatter plot when using `plot = TRUE`. Now correctly adds this to details list when `details = TRUE`.

* `new sampling method` - Added `sample_nc()` based on the algorithm described in [Melville & Stone (2016)](https://doi.org/10.1080/00049158.2016.1218265).

    * This algorithm uses kmeans clustering where the number of clusters is equal to the desired number of samples. Cluster centers are located, which then prompts the nearest neighbour raster pixel for each cluster to be located (assuming default `k` parameter). These nearest neighbours are the output samples. Visualization of the centers and samples can be dispayed if `details = TRUE` is used and `$kplot` is plotted.

* `fixed` - `sample_systematic()` - now has inherent randomness for lower left corner of the tessellation.

* `fixed` - `strat_kmeans()` - solved issue where only first raster layer was being involved in stratification.

# sgsR 1.1.0

* `enhanced` - `sample_strat()` - added parameter `method` that allows users to choose between `"Queinnec"` (default method implemented in previous sgsR versions) and `"random"` (stratified random sampling). The random method ignores much of the functionality of the algorithm to allow users to use standard stratified random sampling approaches without the use of a focal window to locate contiguous stratum cells.

* `fixed` - `sample_strat()` - factor handling improvement - GitHub issue #18.

* `enhanced` -  `calculate_allocation()` - improved documentation for output data frame to make attributes more clear.

* `fixed` - `calculate_representation()` - will now not plot bar chart twice & `NA` values in existing will not be removed.

* `fixed` - `existing` samples with other attributes will now not break sampling using `sample_ahels() / sample_clhs()` if values are `NA`. Variables are also added back to the sample output.

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

* Changes to `terra::distance()` & `terra::classify()` required slight modifications to `calculate_distance()` and `strat_breaks()`.

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
