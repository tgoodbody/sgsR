% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_srs.R
\name{sample_srs}
\alias{sample_srs}
\title{Simple random sampling}
\usage{
sample_srs(
  raster,
  nSamp,
  mindist = NULL,
  access = NULL,
  buff_inner = NULL,
  buff_outer = NULL,
  plot = FALSE,
  filename = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{raster}{spatRaster. Raster to be used for random sampling.}

\item{nSamp}{Numeric. Number of desired samples.}

\item{mindist}{Numeric. Minimum allowable distance between selected
samples. \code{Default = NULL}.}

\item{access}{sf 'LINESTRING' or 'MULTILINESTRING'. Access network.}

\item{buff_inner}{Numeric. Inner buffer boundary specifying distance
from access where plots cannot be sampled.}

\item{buff_outer}{Numeric. Outer buffer boundary specifying distance
from access where plots can be sampled.}

\item{plot}{Logical. Plots output strata raster with samples.}

\item{filename}{Character. Path to write output samples.}

\item{overwrite}{Logical. Choice to overwrite existing \code{filename} if it exists.}
}
\value{
An sf object with \code{nSamp} randomly sampled points.
}
\description{
Randomly sample within a stratification raster extent.
}
\examples{
#--- Load raster and access files ---#
r <- system.file("extdata", "sraster.tif", package = "sgsR")
sr <- terra::rast(r)

a <- system.file("extdata", "access.shp", package = "sgsR")
ac <- sf::st_read(a)

#--- perform simple random sampling ---#
sample_srs(
  raster = sr,
  nSamp = 200,
)

}
\seealso{
Other sample functions: 
\code{\link{sample_ahels}()},
\code{\link{sample_balanced}()},
\code{\link{sample_clhs}()},
\code{\link{sample_existing}()},
\code{\link{sample_nc}()},
\code{\link{sample_strat}()},
\code{\link{sample_sys_strat}()},
\code{\link{sample_systematic}()}
}
\author{
Tristan R.H. Goodbody & Martin Queinnec
}
\concept{sample functions}
