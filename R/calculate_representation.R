#' Compare representation of samples within sraster strata
#'
#' @details Calculate how well sraster strata are represented in existing samples
#'
#' @family calculate functions
#'
#' @inheritParams sample_strat
#' 
#' @return Returns a data.frame of:
#' \itemize{
#' \item{strata} - \code{sraster} strata ID.
#' \item{srasterFreq} - Coverage frequency percent of \code{sraster} strata.
#' \item{sampleFreq} - Sampling frequency percent within \code{sraster} strata.
#' \item{diffFreq} - Difference between \code{srasterFreq} & \code{sampleFreq}. Positive values indicate over representation
#' \item{nSamp} - Number of samples within each strata in \code{existing}.
#' \item{need} - \code{srasterFreq * sum(nSamp)}. Total theoretical number of required samples to be representative of strata coverage.
#' This values is rounded. It is important for the user to consider \code{diffFreq}. A small difference - e.g. 1 percent -
#' in sample vs. sraster frequency could correspond to the algorithm allocating or removing samples that could likely be ignored.
#' }
#'
#' @examples
#' ###--- generate example stratification ---###
#' 
#' #--- load ALS metrics from sgsR internal data ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' 
#' #--- read ALS metrics using the terra package ---#
#' mraster <- terra::rast(r)
#' 
#' #--- perform stratification ---#
#' sraster <- strat_kmeans(mraster = mrasters$zq90,
#'                         nStrata = 6,
#'                         plot = TRUE)
#' 
#' ###--- create existing sample network ---###
#' 
#' #--- load ALS metrics from sgsR internal data ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' 
#' #--- read ALS metrics using the terra package ---#
#' mraster <- terra::rast(r)
#' 
#' #--- simple random sampling ---#
#' existing <- sample_srs(raster = mraster$zq90,
#'                        nSamp = 100)
#'                        
#' #--- calculate representation ---#
#' 
#' calculate_representation(sraster = sraster, 
#'                          existing = existing, 
#'                          plot = TRUE)
#'                          
#' @author Tristan R.H. Goodbody, Martin Queinnec
#'
#' @export

calculate_representation <- function(sraster,
                                     existing,
                                     plot = FALSE) {
  
  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!stringr::str_detect(names(sraster), "strata")) {
    stop("A layer name containing 'strata' does not exist within 'sraster'.")
  }
  
  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object", call. = FALSE)
  }


  ###--- evaluate sample ---###
  #--- determine crs of input sraster ---#
  
  crs <- terra::crs(sraster, proj = TRUE)
  
  #--- extract covariates data from mraster ---#
  
  vals <- terra::as.data.frame(sraster, xy = TRUE, row.names = FALSE) %>%
    dplyr::rename(
      X = x,
      Y = y
    )
  
  #--- Remove NA / NaN / Inf values - calculate frequency of strata coverage ---#
  
  vals_mat <- vals %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(cnt = dplyr::n()) %>%
    dplyr::mutate(srasterFreq = round(cnt / sum(cnt), 2)) %>% 
    dplyr::arrange(desc(srasterFreq))
  
  #--- existing ---#
  existing_mat <- extract_strata(sraster = sraster, existing = existing, data.frame = TRUE) %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(nSamp = dplyr::n()) %>%
    dplyr::mutate(sampleFreq = round(nSamp / sum(nSamp), 2)) %>% 
    dplyr::arrange(desc(sampleFreq))
  
  #--- compare ---#
  rep <- dplyr::left_join(vals_mat,existing_mat, by = "strata") %>%
    replace(is.na(.), 0) %>% # if no samples are within a strata replace NA with 0
    dplyr::mutate(diffFreq = sampleFreq - srasterFreq) %>%
    dplyr::select(strata,srasterFreq,sampleFreq,diffFreq, nSamp) %>%
    dplyr::mutate(need = ceiling(srasterFreq * sum(nSamp)) - nSamp) %>%
    dplyr::arrange(strata)
  
  #--- present barchart if desired ---#
  if (isTRUE(plot)) {
    p <- rep %>%
      dplyr::select(strata, sraster = srasterFreq, samples = sampleFreq) %>%
      tidyr::pivot_longer(c(2,3)) %>%
      ggplot2::ggplot(ggplot2::aes(x = as.factor(strata), y = value, fill = name))+
      ggplot2::geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values=c("#141414", "#5c5c5c"))+
      ggplot2::labs(x = "Strata",
           y = "Frequency",
           title = "Sample representation by strata",
           subtitle = "Strata coverage frequency vs. sampling frequency within strata")+
      theme(legend.position="bottom",
            legend.title=ggplot2::element_blank())
          
    print(p)
  }

return(rep)

}
