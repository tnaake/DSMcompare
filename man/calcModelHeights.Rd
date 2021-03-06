% Generated by roxygen2 (4.0.2): do not edit by hand
\name{calcModelHeights}
\alias{calcModelHeights}
\title{calcModelHeights}
\usage{
calcModelHeights(coordinates, model, method = c("2D", "IDW"),
                         idw = list("p" = 2, "m", "rad" = 5))
}
\arguments{
\item{coordinates}{a matrix or data.frame containing x- and y-coordinates
model a list of \code{data.frame}s with model values with x-coordinates,
y-coordinates, z-coordinates and classes
method character, use "2D" or "IDW"
idw a list with parameters for IDW interpolarion}
}
\value{
a list of data frames with x-coordinates, y-coordinates, z-coordinates
and class.
}
\description{
a function to interpolate z-coordinates at given x- and y-coordinates
}
\details{
see Shepard for IDW.
}
\examples{
\dontrun{
calcModelHeights(coordinates, model, method = c("2D", "IDW"),
                         idw = list("p" = 2, "m", "rad" = 5))
}
}
\author{
Thomas Naake <naake@stud.uni-heidelberg.de>
}

