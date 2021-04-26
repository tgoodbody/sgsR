# sgsR - structurally guided sampling using ALS metrics

Development of sgsR was made possible thanks to the financial support of _______________________.

## Installation

```
install.packages("devtools")
devtools::install_github("https://github.com/tgoodbody/sgsR", build_vignettes = FALSE)
library(sgsR)
```

If you want to build the vignette use `build_vignette = TRUE`.

## Example usage

Stratify ALS metrics using a variety of stratification methods to produce an output stratification raster

```
#--- perform stratification using k-means ---#
kmeans <- strat_kmeans(mraster = mraster, 
                       nstrata = 4)

pcomp <- strat_pcomp(mraster = mraster, 
                      nstrata = 4)

metrics <- strat_metrics(mraster = wall_poly, 
                          metric = "metric", 
                          metric2 = "metric2", 
                          nstrata = 10, 
                          nstrata2 = 5)

```

Sample within stratification rasters

### Simple random sampling

```
srs_wo <- sample_srs(sraster = sraster,
                     n = 50,
                     plot = TRUE)
```

### Stratified sampling

```
#--- stratified sampling with access and buffers provided ---#
strat_w_a <- sample_strat(sraster = sraster,
                          n = 200, 
                          mindist = 200,
                          access = roads,
                          buff_inner = 50,
                          buff_outer = 200,
                          plot = TRUE)
```










