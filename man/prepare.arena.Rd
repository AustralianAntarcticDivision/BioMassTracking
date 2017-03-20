\name{prepare.arena}
\alias{prepare.arena}
\title{Prepares the arena for an estimation of movement}
\description{
This function prepares the arena, i.e. the flow field and the polygon structure,
for use by the other functions in package BiomassTracking.
}
\usage{
prepare.arena(file_U = NULL, file_V = NULL, file_S = NULL, n = dim(file_U)[2], 
              m = dim(file_U)[1], na.string = "-9999", eddy = FALSE,
              smooth = FALSE, noise_sd = 0)
}
\arguments{
  \item{file_U}{Describes the U data, see details.}
  \item{file_V}{Describes the V data, see details.}
  \item{file_S}{Describes the S data, see details.}
  \item{n}{Zonal dimension of the arena, if needed.}
  \item{m}{Meridional dimension of the arena, if needed.}
  \item{na.string}{The character string denoting NA in the input files.}
  \item{eddy}{If TRUE, an eddy is added in the middle of the arena.}
  \item{smooth}{If TRUE, the flow field is smoothed with a 3x3 moving average.}
  \item{noise_sd}{The standard deviation of Gaussian noise added to the flow field.}
}
\details{
The data descriptors \code{file\_U}, \code{file\_V} and \code{file\_S} can come in a number of different 
formats. If they are filenames, the respective data is loaded from these files. In case the data is 
already loaded, they can also be data.frames. If they are NULL, a default construction is invoked, 
where the values \code{n} and \code{m} are needed and the flow is straight in west-east direction. 
They can also be matrices giving the respective data directly.

In the usual case, the data is loaded from files and then prepared to be used. The matrices \code{U} 
and \code{V} are interpolated in every NA position that is not land, while the velocities on land are 
set to 0. The component perpendicular to land of the sea velocities directly at the land-sea boundary 
is set to 0 if it points towards land, so that particles wouldn't cross the boundary if an ideal particle tracking 
algorithm was used.

The data is then rescaled, as the files are supposed to be in m/s format while the particle tracking
algorithms act on a longitude/latitude grid. Therefore, the values of 'U' and 'V' are rescaled to
degree lon/day and degree lat/day, respectively. Note that this rescaling includes a relative increase
of the southern U-velocities compared to the northern U-velocities, because the longitudinal distances
are smaller in the South. Also note that the whole algorithm is meant to work in the southern hemisphere
and has to be adapted to the Northern Hemisphere if that is wanted.
}
\value{
An object of class \code{arena}, describing the arena in which the function is too be used. It is a list containing
elements \code{lat}, \code{lon}, \code{U}, \code{V} and \code{S} (in that order).
  \item{lon, lat}{ Vectors storing the latitude and longitude values of the grid points used.}
  \item{U, V}{Matrices with the corresponding flow velocities in west-east and south-north direction.}
  \item{S}{A matrix in which each grid points has an integer number, either giving the polygon it belongs to 
(if > 0) or stating that this grid point lies on land (if == 0).}
}
\section{Warning}{This function is written only for the Southern Hemisphere, but can be adapted easily to the Northern Hemisphere.}
\references{TODO: my report}
\author{ Thorsten Lenser }


\seealso{\code{\link{particle.tracking}}, \code{\link{biomass.tracking}}, \code{\link{particle.tracking.compare}}, \code{\link{biomass.compare}}, \code{\link{plot.arena}}}
\examples{
data(Udata)
data(Vdata)
data(Sdata)
arena = prepare.arena(Udata,Vdata,Sdata)
plot(arena)
}
\keyword{ misc }% at least one, from doc/KEYWORDS

