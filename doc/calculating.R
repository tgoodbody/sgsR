## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE,results=FALSE-----------------------------
library(sgsR)
library(terra)

#--- Load mraster and access files ---#
r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)

a <- system.file("extdata", "roads.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a, quiet = TRUE)

#--- apply kmeans algorithm to metrics raster ---#
sraster <- strat_kmeans(mraster = mraster, 
                        nStrata = 4) # algorithm will plot output

#--- apply kmeans algorithm to metrics raster ---#
existing <- sample_srs(raster = mraster, # use sraster as input for sampling
                       nSamp = 200, # request 200 samples be taken
                       mindist = 100) # algorithm will plot output


## ----warning=F,message=F------------------------------------------------------
calculate_distance(raster = sraster, # input
                   access = access, # define access road network
                   plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
calculate_pcomp(mraster = mraster, # input
                nComp = 5, # number of components to output
                plot = TRUE, # plot
                details = TRUE) # details about the principal component analysis appended


## ---- warning = F-------------------------------------------------------------
#--- determine sample size based on relative standard error (rse) of 1% ---#
calculate_sampsize(mraster = mraster,
                   rse = 0.01)


## ---- warning = FALSE---------------------------------------------------------
#--- change default threshold sequence values ---# 
#--- if increment and rse are not divisible the closes value will be taken ---#
p <- calculate_sampsize(mraster = mraster,
                   rse = 0.025,
                   start = 0.01,
                   end = 0.08,
                   increment = 0.01,
                   plot = TRUE)

p

## ----warning=F,message=F------------------------------------------------------
#--- perform grid sampling ---#
calculate_allocation(sraster = sraster, 
                     nSamp = 200)

## ----warning=F,message=F------------------------------------------------------
#--- calculate existing samples to include ---#
e.sr <- extract_strata(sraster = sraster, 
                       existing = existing)

calculate_allocation(sraster = sraster, 
                     nSamp = 200, 
                     existing = e.sr)

## ---- warning=F,message=F-----------------------------------------------------
calculate_allocation(sraster = sraster, # stratified raster
                     nSamp = 200, # desired sample number
                     existing = e.sr, #existing samples
                     allocation = "optim", # optimal allocation
                     mraster = mraster$zmax, # metric raster
                     force = TRUE) # force nSamp number


## -----------------------------------------------------------------------------
calculate_allocation(sraster = sraster, # stratified raster
                     nSamp = 20, # desired sample number
                     allocation = "equal") # optimal allocation

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  calculate_coobs(mraster = mraster, # input
#                  existing = existing, # existing samples
#                  cores = 4, # parallel cores to use
#                  details = TRUE, # provide details from algorithm output
#                  plot = TRUE, # plot
#                  filename = tempfile(fileext = ".tif")) # write output raster to tif

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- by default all statistical data are calculated ---#
#  calculate_lhsPop(mraster = mraster) # input

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- statistical analyses can be chosen by setting their parameter to `FALSE` ---#
#  calculate_lhsPop(mraster = mraster, # input
#                   nQuant = 10, # desired number of quantiles
#                   PCA = FALSE) # choose not to calculate PCA's

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- calculate lhsPop details ---#
#  poplhs <- calculate_lhsPop(mraster = mr)
#  
#  calculate_lhsOpt(popLHS = poplhs)

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  calculate_lhsOpt(popLHS = poplhs,
#                   PCA = FALSE,
#                   iter = 200)

