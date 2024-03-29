% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_existing.R
\name{sample_existing}
\alias{sample_existing}
\title{Sample existing}
\usage{
sample_existing(
  existing,
  nSamp,
  raster = NULL,
  type = "clhs",
  access = NULL,
  buff_inner = NULL,
  buff_outer = NULL,
  details = FALSE,
  filename = NULL,
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{existing}{sf 'POINT'.  Existing plot network.}

\item{nSamp}{Numeric. Number of desired samples.}

\item{raster}{SpatRaster. Raster to guide the location of the samples. If \code{type = "clhs"} this raster can also
be used to define the population distributions to be used for sampling.}

\item{type}{Character. A string indicating the type of sampling method to use.
Possible values are \code{"clhs"}, \code{"balanced"}, \code{"srs"} and \code{"strat"}.}

\item{access}{sf. Road access network - must be lines.}

\item{buff_inner}{Numeric. Inner buffer boundary specifying distance
from access where plots cannot be sampled.}

\item{buff_outer}{Numeric. Outer buffer boundary specifying distance
from access where plots can be sampled.}

\item{details}{Logical. If \code{FALSE} (default) output is sf object of
systematic samples. If \code{TRUE} returns a list of sf objects where \code{tessellation}
is the tessellation grid for sampling, and \code{samples} are the systematic samples.}

\item{filename}{Character. Path to write output samples.}

\item{overwrite}{Logical. Choice to overwrite existing \code{filename} if it exists.}

\item{...}{Additional arguments for the sampling method selected.}
}
\value{
An sf object of samples or a list object if \code{details = TRUE}
}
\description{
Sub-sample an existing sample. Four sampling methods are available:
\code{clhs}, \code{balanced}, \code{srs} and \code{strat}.
}
\note{
When \code{type = "clhs"} or \code{type = "balanced"} all attributes in \code{existing} will be used for sampling.
Remove attributes not indented for sampling' prior to using this algorithm.
}
\examples{
r <- system.file("extdata", "mraster.tif", package = "sgsR")
mr <- terra::rast(r)

#--- generate an existing sample adn extract metrics ---#
e <- sample_systematic(raster = mr, cellsize = 200)
e <- extract_metrics(existing = e, mraster = mr)

#--- perform clhs (default) sub-sampling ---#
sample_existing(
  existing = e,
  nSamp = 50
)

#--- perform balanced sub-sampling ---#
sample_existing(
  existing = e,
  nSamp = 50,
  type = "balanced"
)

#--- perform simple random sub-sampling ---#
sample_existing(
  existing = e,
  nSamp = 50,
  type = "srs"
)

}
\seealso{
Other sample functions: 
\code{\link{sample_ahels}()},
\code{\link{sample_balanced}()},
\code{\link{sample_clhs}()},
\code{\link{sample_nc}()},
\code{\link{sample_srs}()},
\code{\link{sample_strat}()},
\code{\link{sample_sys_strat}()},
\code{\link{sample_systematic}()}
}
\author{
Tristan R.H. Goodbody
}
\concept{sample functions}
