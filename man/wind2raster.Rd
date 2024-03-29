% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wind_functions2.R
\name{wind2raster}
\alias{wind2raster}
\title{Wind-data to raster file}
\usage{
wind2raster(x)
}
\arguments{
\item{x}{an "rWind list" obtained by wind.fit}
}
\value{
A raster stack or a list of raster stacks representing wind direction
and speed.
}
\description{
wind2raster crates a raster stack (gridded) with 2 layers: wind speed and
wind direction for an object of \code{rWind}.
Latitude and longitude values are used to locate raster file and to create
raster using rasterFromXYZ function from raster package. If the input file is
a list of wind data created by wind.dl, a list of raster stacks will be
returned
}
\details{
WGS84 datum (non-projected) CRS is selected by default to build the raster
file.
}
\examples{

data(wind.data)

# Create raster stack from the downloaded data with wind directon and speed
# layers

wind2raster(wind.data)
}
\seealso{
\code{\link{wind.dl}}
}
\author{
Javier Fernández-López (jflopez.bio@gmail.com)
}
\keyword{~gfs}
\keyword{~wind}
