% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strat_quantiles.R
\name{strat_quantiles}
\alias{strat_quantiles}
\title{Quantiles stratification}
\usage{
strat_quantiles(
  mraster,
  nStrata,
  map = FALSE,
  plot = FALSE,
  details = FALSE,
  filename = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{mraster}{spatRaster. ALS metrics raster.}

\item{nStrata}{Numeric. A positive integer representing the number of desired classes
or a numeric vector of probabilities with values between \code{0-1}. If \code{mraster} has multiple layers,
\code{nStrata} must be a list with an equal number of objects.}

\item{map}{Logical. Map individual stratified layers to a combined stratification. Will output a multi-layer
\code{SpatRaster} with individual stratifications for each \code{mraster} layer and an additional mapped stratification
 named \code{"strata"}.}

\item{plot}{Logical. Plots output strata raster and visualized
strata with boundary dividers.}

\item{details}{Logical. If \code{FALSE} (default) output is only
stratification raster. If \code{TRUE} return a list where \code{$details} is additional
stratification information and \code{$raster} is the output stratification spatRaster.}

\item{filename}{Character. Path to write stratified raster to disc.}

\item{overwrite}{Logical. Specify whether \code{filename} should be overwritten on disc.}
}
\value{
Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.

When a list is returned:
\enumerate{
\item \code{details} lookUp table for stratification(s).
\item \code{raster} is a stratified \code{spatRaster} based on quantiles
\item \code{plot} is a \code{ggplot} histogram / scatter plot object (depends on whether metric2 was supplied).
Histogram shows distribution and break points.
}
}
\description{
Stratify metric raster using metric quantiles.
}
\examples{
#--- Load raster and existing plots---#
r <- system.file("extdata", "mraster.tif", package = "sgsR")
mr <- terra::rast(r)

strat_quantiles(
  mraster = mr$zq90,
  nStrata = 4
)

strat_quantiles(
  mraster = mr[[1:2]],
  nStrata = list(c(0.2, 0.4, 0.8), 3),
  map = TRUE
)

}
\seealso{
Other stratify functions: 
\code{\link{strat_breaks}()},
\code{\link{strat_kmeans}()},
\code{\link{strat_map}()},
\code{\link{strat_poly}()}
}
\author{
Tristan R.H. Goodbody
}
\concept{stratify functions}
