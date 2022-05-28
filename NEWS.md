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
