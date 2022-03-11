## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F------------------------------------------------------
library(sgsR)
library(terra)
library(sf)

#--- Load mraster and access files ---#
r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)

## ----warning=F,message=F------------------------------------------------------
a <- system.file("extdata", "roads.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a)

## ----warning=F,message=F------------------------------------------------------
terra::plot(mraster[[1]])
terra::plot(access, add = TRUE, col = "black")


## ----warning=F,message=F------------------------------------------------------
#--- apply kmeans algorithm to metrics raster ---#
sraster <- strat_kmeans(mraster = mraster, # use mraster as input for sampling
                        nStrata = 4, # algorithm will produce 4 strata
                        plot = TRUE) # algorithm will plot output


## ----warning=F,message=F------------------------------------------------------
#--- set seed ---#
set.seed(2021)

#--- apply kmeans algorithm to metrics raster ---#
existing <- sample_srs(raster = mraster, # use mraster as input for sampling
                       nSamp = 200, # request 200 samples be taken
                       mindist = 100, # define that samples must be 100 m apart
                       plot = TRUE) # algorithm will plot output


## ----pipe, eval= FALSE--------------------------------------------------------
#  #--- non piped ---#
#  sraster <- strat_kmeans(mraster = mraster, # use mraster as input for sampling
#                          nStrata = 4, # algorithm will produce 4 strata
#                          plot = TRUE) # algorithm will plot output
#  
#  existing <- sample_srs(raster = sraster, # use mraster as input for sampling
#                         nSamp = 200, # request 200 samples be taken
#                         mindist = 100, # define that samples must be 100 m apart
#                         plot = TRUE) # algorithm will plot output
#  
#  extract_metrics(mraster = mraster,
#                  existing = existing)
#  
#  
#  #--- piped ---#
#  strat_kmeans(mraster = mraster, nStrata = 4) %>%
#    sample_srs(., nSamp = 200, mindist = 100) %>%
#    extract_metrics(mraster = mraster, existing = .)
#  

