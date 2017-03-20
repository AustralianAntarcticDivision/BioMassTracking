## Visualise a flow field as a field of arrows.

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