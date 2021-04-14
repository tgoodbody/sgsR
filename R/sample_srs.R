sample_srs <- function(raster,
                       n){
  
  crs <- crs(raster)
  
  cells <- terra::spatSample(raster, n, "random", cells = TRUE)
  
  v <- raster[cells]
  
  xy <- terra::xyFromCell(raster,cells)
  
  samples <- cbind(xy,v)
  
  coords <- xy %>% 
    as.data.frame() %>% 
    st_as_sf(., coords = c("x","y"))
  
  st_crs(coords) <- crs
  
  plot(raster[[1]])
  plot(coords,add=T)
  
  return(samples)
  
  
}
