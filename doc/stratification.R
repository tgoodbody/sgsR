## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE,results=FALSE-----------------------------
library(sgsR)
library(terra)
library(sf)

#--- Load mraster and access files ---#
r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)

a <- system.file("extdata", "roads.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a, quiet = TRUE)

#--- apply kmeans algorithm to metrics raster ---#
sraster <- strat_kmeans(mraster = mraster, # use mraster for stratification
                        nStrata = 4) # algorithm will plot output

#--- apply kmeans algorithm to metrics raster ---#
existing <- sample_srs(raster = sraster, # use sraster as input for sampling
                       nSamp = 200, # request 200 samples be taken
                       mindist = 100) # algorithm will plot output


## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using k-means ---#
strat_kmeans(mraster = mraster, # input
             nStrata = 5) # algorithm will produce 4 strata


## ----warning=F,message=F------------------------------------------------------
strat_kmeans(mraster = mraster, # input
             nStrata = 10, # algorithm will produce 10 strata
             iter = 1000, # set minimum number of iterations to determine kmeans centers
             algorithm = "MacQueen", # use MacQueen algorithm
             plot = TRUE) # plot output

## ----warning=F,message=F------------------------------------------------------
strat_kmeans(mraster = mraster, # input
             nStrata = 5, # algorithm will produce 4 strata
             center = FALSE, # do not center data
             scale = FALSE, # do not scale data
             plot = TRUE, # plot output
             filename = tempfile(fileext = ".tif"), # write output sraster to file
             overwrite = TRUE) # overwrite file on disc if it exists


## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using principal components ---#
strat_pcomp(mraster = mraster, # input
            nStrata = 5, # 5 strata with primary PC only
            plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
strat_pcomp(mraster = mraster, # input
            nStrata = 4, # 4 strata with primary
            nStrata2 = 4, # 4 strata with secondary PC - will produce 16 output strata
            plot = TRUE) # produce output details

## ----warning=F,message=F------------------------------------------------------
strat_pcomp(mraster = mraster, # input
            nStrata = 3, # 3 strata with primary PC
            nStrata2 = 3, # 4 strata with secondary PC - will produce 9 output strata
            filename = tempfile(fileext = ".tif")) # write output sraster to file

## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using user-defined breaks ---#

#--- define breaks for metric ---#
breaks <- c(seq(0,100,20))

breaks

#--- perform stratification using user-defined breaks ---#

values <- terra::values(mraster$zmax)

#--- define breaks for metric ---#
breaks2 <- quantile(values, na.rm=TRUE)

breaks2


## ----warning=F,message=F------------------------------------------------------
#--- stratify on 1 metric only ---#

strat_breaks(mraster = mraster$zmean,
             breaks = breaks,
             details = TRUE,
             plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
#--- stratify on 1 metric only ---#

strat_breaks(mraster = mraster$zmax,
             breaks = breaks2,
             details = TRUE,
             plot = TRUE)

## -----------------------------------------------------------------------------
#--- load in polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")

fri <- sf::st_read(poly)

## -----------------------------------------------------------------------------
#--- stratify polygon coverage ---#
#--- specify polygon attribute to stratify ---#

attribute <- "NUTRIENTS"

#--- specify features within attribute & how they should be grouped ---#
#--- as a single vector ---#

features <- c("poor", "rich", "medium")

srasterpoly <- strat_poly(poly = fri, # input polygon
                          attribute = attribute, # attribute to stratify by
                          features = features, # features within attribute
                          raster = sraster, # raster to define extent and resolution for output
                          plot = TRUE) # plot output

## -----------------------------------------------------------------------------
#--- or as multiple lists ---#
g1 <- "poor"
g2 <- c("rich", "medium")

features <- list(g1, g2)

strat_poly(poly = fri,
           attribute = attribute,
           features = features,
           raster = sraster,
           plot = TRUE,
           details = TRUE)

## -----------------------------------------------------------------------------
#--- map srasters ---#
strat_map(sraster = srasterpoly, # strat_poly 3 class stratification
          sraster2 = sraster, # strat_kmeans 4 class stratification
          plot = TRUE)


## -----------------------------------------------------------------------------
strat_map(sraster = srasterpoly, # strat_poly 3 class stratification
          sraster2 = sraster, # strat_poly 3 class stratification
          stack = TRUE, # stack input and oputput strata into multi layer ouput raster
          details = TRUE, # provide additional details
          plot = TRUE) # plot output

