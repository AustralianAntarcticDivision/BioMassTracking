## Visualise a flow field as a field of arrows.



#' Plot a flow field
#' 
#' This function plots a flow field given as west-east and south-north velocity
#' values at regularly spaced grid points.
#' 
#' 
#' @param u The west-east velocities as a matrix.
#' @param v The south-north velocities as a matrix.
#' @param xpos Positions of the arrowbases in x-direction.
#' @param ypos Positions of the arrowbases in y-direction.
#' @param add If TRUE, the flow field is added to the current plot.
#' @param scale The length of the longest arrow is 'scale' * 'distance between
#' grid points'.
#' @author Thorsten Lenser
#' @seealso \code{\link{plot.arena}}
#' @keywords misc
#' @examples
#' 
#' data(Udata,Vdata,Sdata)
#' arena = prepare.arena(Udata,Vdata,Sdata)
#' 
#' # First use quiver only
#' quiver(arena$U,arena$V)
#' 
#' # Now use larger arrows, looks better in this case
#' quiver(arena$U,arena$V,scale=5)
#' 
#' # Together with the polygon structure
#' plot(arena,arrow_scale=5)
#' 
#' 
quiver <- function(u,v,xpos=col(u),ypos=row(u),add=FALSE,scale=1)
{
    tmp_u <- par("usr")
    tmp_p <- par("pin")
    units = c(tmp_p[1]/(tmp_u[2] - tmp_u[1]), tmp_p[2]/(tmp_u[4] - tmp_u[3]))

    if (add == TRUE){ 
      matplot(xpos,ypos,type="p",cex=0,add=TRUE)
    }
    else {
      matplot(xpos,ypos,type="p",cex=0)
    }
	
    unit = min(units)

    speed <- sqrt(u*u+v*v)
    maxspeed <- max(speed,na.rm=T)

    u <- u*scale/maxspeed
    v <- v*scale/maxspeed

    # No warnings needed here, they just relate to zero-length arrows
    # that are skipped in the drawing.
    suppressWarnings(arrows(xpos,ypos,xpos+u,ypos+v,length=unit/5,code=2))
} 
