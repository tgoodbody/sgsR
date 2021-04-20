### sgsR vignette ###

#--- import ALS metrics raster ---#

raster <- rast("C:/Users/goodb/Documents/UBC/post_doc/sgsR/vig/data/wall_metrics_small.tif")

#--- import forest inventory polygon and mask unwanted areas ---#
poly <- vect("C:/Users/goodb/Documents/UBC/post_doc/sgsR/vig/data/inventory_polygons.shp")

poly_subset <- poly[poly$POLYTYPE == "FOR" & poly$OWNER == 1, ]
poly_subset <- terra::aggregate(poly_subset,dissolve = TRUE)

#--- mask input ALS raster using polygon layer ---#

wall_poly <- terra::mask(raster,poly_subset)

#--- import access layer to be used during sampling if desired ---#

roads <- vect("C:/Users/goodb/Documents/UBC/post_doc/sgsR/vig/data/roads.shp")

#--- perform stratification using k-means ---#
kmeans <- strat_kmeans(raster = wall_poly, k = 8) ### note some values k dont seem to give the right output

#--- perform stratification using OSB ---#
#--- note that this one can take a while ---#
osb <- strat_osb(raster = wall_poly, metric = "wal_5", h = 4, n = 100) # should integrate functionality to do 2 metrics concurrently and create unique strata

#--- perform stratification using principal components ---#
pcomp <- strat_pcomp(raster = wall_poly, ncp = 2, b1 = 4, b2 = 3) # should integrate functionality to allow users to define strata based on more than just PC1 and PC2

#--- perform stratification using individual metrics ---#
metrics <- strat_metrics(raster = wall_poly, metric = "wal_5", metric2 = "wal_2", b = 10, b2 = 5)

#--- define desired stratification raster ---#
raster <- osb$raster

#--- sampling **without** access defined---#

srs_wo <- sample_srs(raster = raster,
                  ns = 50,
                  mindist = 200)

strat_wo <- sample_strat(raster = raster,
                      ns = 50, 
                      mindist = 200)

#--- sampling **with** access defined---#

srs_w <- sample_srs(raster = raster,
                  ns = 200,
                  mindist = 200,
                  access = roads,
                  buff_inner = 50,
                  buff_outer = 200)

strat_w <- sample_strat(raster = raster,
                      ns = 200, 
                      mindist = 200, 
                      access = roads, 
                      buff_inner = 50,
                      buff_outer = 200,
                      buff_extend = 100,
                      buff_max = 600)

#--- extract strata from raster for already existing sample network ---#
#--- we use random samples defined above ---#

existing <- extract_existing(raster,srs_wo)

#--- sampling **with** access defined **and** existing samples defined ---#

strat_w_e <- sample_strat(raster = raster,
                        ns = 200, 
                        mindist = 200, 
                        access = roads, 
                        existing = existing,
                        buff_inner = 50,
                        buff_outer = 200,
                        buff_extend = 100,
                        buff_max = 600)

#--- extract metrics from multi-band ALS raster for potential modeling ---#

metrics <- extract_metrics(wall_poly,strat_w_e$samples)


