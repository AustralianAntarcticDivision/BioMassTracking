\name{quiver}
\alias{quiver}
\title{Plot a flow field}
\description{
This function plots a flow field given as west-east and south-north velocity values
at regularly spaced grid points. 
}
\usage{
quiver(u, v, xpos = col(u), ypos = row(u), add = FALSE, scale=1)
}
\arguments{
  \item{u}{The west-east velocities as a matrix.}
  \item{v}{The south-north velocities as a matrix.}
  \item{xpos}{Positions of the arrowbases in x-direction.}
  \item{ypos}{Positions of the arrowbases in y-direction.}
  \item{add}{If TRUE, the flow field is added to the current plot.}
  \item{scale}{The length of the longest arrow is 'scale' * 'distance between grid points'.}
}

\author{ Thorsten Lenser }

\seealso{\code{\link{plot.arena}}}
\examples{
data(Udata,Vdata,Sdata)
arena = prepare.arena(Udata,Vdata,Sdata)

# First use quiver only
quiver(arena$U,arena$V)

# Now use larger arrows, looks better in this case
quiver(arena$U,arena$V,scale=5)

# Together with the polygon structure
plot(arena,arrow_scale=5)

}
\keyword{ misc }% at least one, from doc/KEYWORDS

