### sgsR vignette ###

#--- import ALS metrics raster ---#

raster <- rast("C:/Users/goodb/Documents/UBC/post_doc/RMFinventory/inst/extdata/wall_metrics_small.tif")

#--- import forest inventory polygon and mask unwanted areas ---#
poly <- vect(system.file("extdata/inventory_polygons","inventory_polygons.shp", package = "RMFinventory"))

poly_subset <- poly[poly$POLYTYPE == "FOR" & poly$OWNER == 1, ]
poly_subset <- terra::aggregate(poly_subset,dissolve = TRUE)

#--- mask input ALS raster using polygon layer ---#

raster <- terra::mask(raster,poly_subset)

#--- import access layer to be used during sampling if desired ---#

roads <- vect("C:/Users/goodb/Documents/UBC/post_doc/RMFinventory/inst/extdata/roads/roads.shp")


#--- perform stratification using k-means ---#
kmeans <- strat_kmeans(raster = raster[[5]], k = 4)

#--- perform stratification using OSB ---#
osb <- strat_osb(raster = raster, metric = "wal_5", h = 4, n = 10)

#--- perform stratification using principal components ---#
pcomp <- strat_pcomp(raster = raster, ncp = 2, b1 = 4, b2 = 3)

#--- perform stratification using individual metrics ---#
metrics <- strat_metrics(raster = raster, metric = "wal_5", metric2 = "wal_2", b = 10, b2 = 5)


#--- sampling without access defined---#

srs_wo <- sample_srs(raster = kmeans$raster,
                  ns = 15,
                  mindist = 200)

strat_wo <- sample_strat(raster = kmeans$raster,
                      ns = 200, 
                      mindist = 200)

#--- sampling with access defined---#

srs_w <- sample_srs(raster = raster,
                  ns = 200,
                  mindist = 200,
                  access = roads,
                  buff_inner = 50,
                  buff_outer = 200)

strat_w <- sample_strat(raster = pcomp$raster,
                      ns = 200, 
                      mindist = 200, 
                      access = roads, 
                      buff_inner = 50,
                      buff_outer = 200,
                      buff_extend = 100,
                      buff_max = 600)

#--- extract strata from raster for already existing sample network ---#

existing <- extract_existing(kmeans$raster,srs_wo)

#--- sampling with access defined and existing samples ---#

strat_w_e <- sample_strat(raster = kmeans$raster,
                        ns = 200, 
                        mindist = 200, 
                        access = roads, 
                        existing = existing,
                        buff_inner = 50,
                        buff_outer = 200,
                        buff_extend = 100,
                        buff_max = 600)


