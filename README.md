
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sgsR - structurally guided sampling <img src="man/figures/logo.png" align="right" width="200" />

<!-- badges: start -->
<!-- badges: end -->

`sgsR` is designed to implement structurally guided sampling approaches
for enhanced forest inventories. The package was designed to function
using rasterized airborne laser scanning (ALS; Lidar) metrics to allow
for stratification of forested areas based on structure.

## Installation

You can install the released version of sgsR from
[Github](https://github.com/tgoodbody/sgsR) with:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/tgoodbody/sgsR")
library(sgsR)
```

## Implementation

-   Describe package fundamentals - `vignette("sgsR")`

-   Overview of sampling algorithms - `vignette("sampling")`

-   Overview of stratification algorithms - `vignette("stratification")`

-   Overview of calculate algorithms - `vignette("calculating")`

## Collaborators

We are thankful for continued collaboration with academic, private
industry, and government institutions to help improve `sgsR`. Special
thanks to to:

| Collaborator                                                                                                  | Affiliation                                                             |
|:--------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------|
| [Martin Queinnec](https://www.researchgate.net/profile/Martin-Queinnec)                                       | University of British Columbia                                          |
| [Joanne C. White](https://scholar.google.ca/citations?user=bqjk4skAAAAJ&hl=en)                                | Canadian Forest Service                                                 |
| [Piotr Tompalski](https://scholar.google.ca/citations?user=RtYdz0cAAAAJ&hl=en)                                | Canadian Forest Service                                                 |
| [Andrew T. Hudak](https://scholar.google.ca/citations?hl=en&user=bdn7YVoAAAAJ)                                | United States Forest Service                                            |
| [Ruben Valbuena](https://scholar.google.com/citations?user=Nx336TQAAAAJ&hl=en)                                | Swedish University of Agricultural Sciences                             |
| [Antoine LeBoeuf](https://scholar.google.com/citations?user=wGsKOK8AAAAJ&hl=en)                               | Ministère des Forêts, de la Faune et des Parcs                          |
| [Ian Sinclair](http://www.infogo.gov.on.ca/infogo/home.html#empProfile/332620/en)                             | Ministry of Northern Development, Mines, Natural Resources and Forestry |
| [Grant McCartney](https://www.signalhire.com/profiles/grant-mccartney%27s-email/99719223)                     | Forsite Consulting                                                      |
| [Jean-Francois Prieur](https://www.researchgate.net/scientific-contributions/Jean-Francois-Prieur-2142960944) | Université du Québec à Montréal Alumni                                  |
| [Murray Woods](https://www.researchgate.net/profile/Murray-Woods)                                             | Ontario Ministry of Natural Resources                                   |

## Funding

Development of sgsR was made possible thanks to the financial support of
the [Canadian Wood Fibre Centre’s Forest Innovation
Program](https://www.nrcan.gc.ca/science-and-data/funding-partnerships/funding-opportunities/forest-sector-funding-programs/forest-innovation-program/13137).
