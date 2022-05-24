
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sgsR - structurally guided sampling <img src="man/figures/logo.png" align="right" width="200" />

<!-- badges: start -->

![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg)
[![R-CMD-check](https://github.com/tgoodbody/sgsR/workflows/R-CMD-check/badge.svg)](https://github.com/tgoodbody/sgsR/actions)
<!-- badges: end -->

`sgsR` is designed to implement structurally guided sampling approaches
for enhanced forest inventories. The package was designed to function
using rasterized airborne laser scanning (ALS; Lidar) metrics to allow
for stratification of forested areas based on structure.

## Installation :computer:

You can install the released version of sgsR from
[Github](https://github.com/tgoodbody/sgsR) with:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/tgoodbody/sgsR")
library(sgsR)
```

## Example usage :bar_chart:

``` r
#--- Load mraster files ---#
r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)

#--- apply kmeans algorithm to mraster ---#
sraster <- strat_kmeans(mraster = mraster, # use mraster as input for stratification
                        nStrata = 4, # produce 4 strata
                        plot = TRUE) # plot output
                        
#--- apply stratified sampling ---#
existing <- sample_strat(sraster = sraster, # use sraster as input for sampling
                         nSamp = 200, # request 200 samples
                         mindist = 100, # samples must be 100 m apart
                         plot = TRUE) # plot output
```

## Vignettes :books:

Check out [the online
documentation](https://tgoodbody.github.io/sgsR/index.html) to see how
`sgsR` functions and how you might use it for your work!

Vignettes include:

-   Package fundamentals - `vignette("sgsR")`

-   Sampling algorithms - `vignette("sampling")`

-   Stratification algorithms - `vignette("stratification")`

-   Calculate algorithms - `vignette("calculating")`

## Collaborators :woman: :man:

We are thankful for continued collaboration with academic, private
industry, and government institutions to help improve `sgsR`. Special
thanks to to:

| Collaborator                                                                                                  | Affiliation                                                                       |
|:--------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
| [Martin Queinnec](https://www.researchgate.net/profile/Martin-Queinnec)                                       | University of British Columbia                                                    |
| [Joanne C. White](https://scholar.google.ca/citations?user=bqjk4skAAAAJ&hl=en)                                | Canadian Forest Service                                                           |
| [Piotr Tompalski](https://scholar.google.ca/citations?user=RtYdz0cAAAAJ&hl=en)                                | Canadian Forest Service                                                           |
| [Andrew T. Hudak](https://scholar.google.ca/citations?hl=en&user=bdn7YVoAAAAJ)                                | United States Forest Service                                                      |
| [Ruben Valbuena](https://scholar.google.com/citations?user=Nx336TQAAAAJ&hl=en)                                | Swedish University of Agricultural Sciences                                       |
| [Antoine LeBoeuf](https://scholar.google.com/citations?user=wGsKOK8AAAAJ&hl=en)                               | Ministère des Forêts, de la Faune et des Parcs                                    |
| [Ian Sinclair](http://www.infogo.gov.on.ca/infogo/home.html#empProfile/332620/en)                             | Ministry of Northern Development, Mines, Natural Resources and Forestry           |
| [Grant McCartney](https://www.signalhire.com/profiles/grant-mccartney%27s-email/99719223)                     | Forsite Consultants Ltd.                                                          |
| [Jean-Francois Prieur](https://www.researchgate.net/scientific-contributions/Jean-Francois-Prieur-2142960944) | Université de Sherbrooke                                                          |
| [Murray Woods](https://www.researchgate.net/profile/Murray-Woods)                                             | (Retired) Ministry of Northern Development, Mines, Natural Resources and Forestry |

## Funding :raised_hands:

Development of sgsR was made possible thanks to the financial support of
the [Canadian Wood Fibre Centre’s Forest Innovation
Program](https://www.nrcan.gc.ca/science-and-data/funding-partnerships/funding-opportunities/forest-sector-funding-programs/forest-innovation-program/13137).
