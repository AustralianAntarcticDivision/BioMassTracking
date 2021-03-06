## Show the arena 
## This is a method for the class 'arena' to extend the generic 
## function 'plot', so a call to 'plot' with an argument
## of class 'arena' will invoke this method.



#' Show the arena
#' 
#' Plots a graphical representation of the flow field together with the polygon
#' structure. This is a method for the generic function \code{plot} and can
#' therefore be called as \code{plot(arena_object)}.
#' 
#' An \code{arena} object describes the arena in which the function is to be
#' used. It is a list containing elements \code{lat}, \code{lon}, \code{U},
#' \code{V} and \code{S} (in that order), where \code{lat} and \code{lon} are
#' vectors storing the latitude and longitude values of the grid points used,
#' \code{U} and \code{V} are matrices with the corresponding flow velocities in
#' west-east and south-north direction, and \code{S} is a matrix in which each
#' grid points has an integer number, either giving the polygon it belongs to
#' (if > 0) or stating that this grid point lies on land (if == 0).
#' 
#' @param x An object of class \code{arena}, see details.
#' @param arrow_scale The length of the longest arrow is \code{arrow\_scale}
#' times the distance between two neighbouring grid points.
#' @param grid If TRUE, a grid is drawn marking the boxes around the grid
#' points.
#' @param \dots Further arguments to the \code{plot} funtion.
#' @author Thorsten Lenser
#' @seealso \code{\link{plot}}, \code{\link{prepare.arena}}
#' @references TODO: my report
#' @keywords misc
#' @examples
#' 
#' data(Udata)
#' data(Vdata)
#' data(Sdata)
#' 
#' arena = prepare.arena(Udata,Vdata,Sdata)
#' plot(arena)
#' @export
plot.arena <- function(x,arrow_scale=1,grid=FALSE,...){

if(class(x) != "arena"){
	stop("First argument has to be an object of class 'arena'!")
}
U = x$U
V = x$V
lat = x$lat
lon = x$lon
S = x$S
m = dim(U)[1]
n = dim(U)[2]

# Plot polygon structure and flow field
image(x=lon,y=lat,t(S),col=topo.colors(12))
suppressWarnings(quiver(U,V,matrix(rep(lon,m),nrow=m,byrow=T),matrix(rep(lat,n),nrow=m,byrow=F),add=TRUE,scale=arrow_scale)) # No warnings, as lots of arrows have length 0
if(grid == TRUE) grid(nx=n,ny=m)

}
