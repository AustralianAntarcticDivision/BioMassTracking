## Show the arena 
## This is a method for the class 'arena' to extend the generic 
## function 'plot', so a call to 'plot' with an argument
## of class 'arena' will invoke this method.

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
