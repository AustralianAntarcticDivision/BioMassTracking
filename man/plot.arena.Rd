\name{plot.arena}
\alias{plot.arena}
\title{Show the arena}
\description{
Plots a graphical representation of the flow field together with 
the polygon structure. This is a method for the generic function
\code{plot} and can therefore be called as \code{plot(arena_object)}.
}
\usage{
plot.arena(x, arrow_scale = 1, grid = FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{arena}, see details.}
  \item{arrow_scale}{The length of the longest arrow is \code{arrow\_scale} times the distance between two neighbouring grid points.}
  \item{grid}{If TRUE, a grid is drawn marking the boxes around the grid points.}
  \item{\dots}{Further arguments to the \code{plot} funtion.}
}
\details{
An \code{arena} object describes the arena in which the function is to be used. It is a list
containing elements \code{lat}, \code{lon}, \code{U}, \code{V} and \code{S} (in that order), where
\code{lat} and \code{lon} are vectors storing the latitude and longitude values of the grid points
used, \code{U} and \code{V} are matrices with the corresponding flow velocities in west-east
and south-north direction, and \code{S} is a matrix in which each grid points has
an integer number, either giving the polygon it belongs to (if > 0) or stating that this grid point
lies on land (if == 0).
}
\references{TODO: my report}
\author{ Thorsten Lenser }

\seealso{\code{\link{plot}}, \code{\link{prepare.arena}}}
\examples{
data(Udata)
data(Vdata)
data(Sdata)

arena = prepare.arena(Udata,Vdata,Sdata)
plot(arena)
}
\keyword{ misc }% at least one, from doc/KEYWORDS

