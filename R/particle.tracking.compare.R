## This version of the particle tracking algorithms is not used to estimate
## parameters, but to compare against a coarser model of biomass movement.
##
## If set to silent, no text or graphical output will appear.



#' Track the movement of biomass as passive drifters
#' 
#' This functions tracks the movement of biomass through the given polygons by
#' following the proportionate distribution of a large number of passive
#' drifter (small fluid parcels) trajectories, which are calculated using a
#' 4th-order Runge-Kutta integration scheme.
#' 
#' There are two different starting conditions for the flow of biomass: A
#' starting distribution can be given or the biomass can flow into the area
#' through a specified inflow polygon in a certain number of time steps.
#' 
#' A call to the C-function \code{particle\_tracking\_compare\_inflow} does the
#' bulk of work, which is a lot faster than the same implementation in R due to
#' the greatly increased efficiency of 'for' loops. It is important to note
#' that in the present form, particles can get caught in 'voids' in the flow
#' field, i.e. in places where their velocity is 0. If this happens for too
#' many particles (a warnings message is given at the end, indicating the
#' number of particles left in the domain), the estimation should be repeated
#' with a higher values of diffusion.
#' 
#' There is considerable space for improvement in replacing the fixed time step
#' Runge-Kutta scheme with a more sophisticated numerical integrator. However,
#' this requires some fiddling with the underlying C functions and has not been
#' tried yet.
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
#' The default value for \code{subdiv} is choosen in such a way that a particle
#' will on average move for 1/10 of the distance between two grid points in
#' each time step. Although a variable time step would be preferable, this
#' yields an accurate estimation in most cases.
#' 
#' Diffusion is included by adding the product of a normal random variable and
#' the mean velocity to the calculated velocity in either direction. The
#' diffusion coefficient gives the standard variation of the normal random
#' number. Diffusion is needed to get particles out of 'voids', i.e. places in
#' the flow field where the adjacent velocities cancel out and the particles
#' get stuck.
#' 
#' @param arena An object of class \code{arena}, see details.
#' @param num_particles The number of particles used.
#' @param t_step The time step in which the results should be presented.
#' @param end_t_counter The number of time steps that should be recorded.
#' @param start_setup A vector giving a starting biomass distribution. Specify
#' either this or the next two arguments.
#' @param infl_poly The polygon into which the biomass enters. Has to have
#' positive inflow from the outside.
#' @param infl_time The number of time steps over which the biomass is to be
#' inserted.
#' @param subdiv The number of Runge-Kutta time steps per recorded time step.
#' See details for the default value.
#' @param diffusion The coefficient of diffusion, see details.
#' @param silent If silent is FALSE, the movement of the biomass will be shown
#' graphically, with the user having to hit 'enter' after every time step
#' @param graphics If graphics is TRUE, the remaining particles at the end will
#' be shown. A plot of of the arena has to be done before.
#' @return A matrix of dimension 'number of polygons' x 'number of time steps',
#' where the last polygon is the 'outside'. Each column holds the biomass
#' distribution in the respective time step.
#' @author Thorsten Lenser
#' @seealso \code{\link{biomass.tracking}},
#' \code{\link{particle.tracking.compare}}, \code{\link{biomass.compare}},
#' \code{\link{prepare.arena}}
#' @references TODO: My report
#' @keywords misc
#' @examples
#' 
#' data(Udata)
#' data(Vdata)
#' data(Sdata)
#' 
#' arena = prepare.arena(Udata,Vdata,Sdata)
#' plot(arena)
#' 
#' # Use more particles if realistic results are needed
#' mk = particle.tracking(arena,400,5000,diffusion=0.6,graphics=TRUE)
#' 
#' # Now estimate the biomass movement
#' N1 = biomass.tracking(mk,seq(0,4900,by=100),infl_poly=2)
#' 
#' # Get a particle tracking result to compare the above to
#' N2 = particle.tracking.compare(arena,400,100,50,diffusion=0.6,infl_poly=2)
#' 
#' # Compare the results
#' biomass.compare(N1,N2)
#' 
particle.tracking.compare <- function(arena,num_particles,t_step,end_t_counter,start_setup=NULL,infl_poly=1,infl_time=10,subdiv=NULL,diffusion=0,silent=TRUE,graphics=FALSE){

if(class(arena) != "arena"){
	stop("First argument has to be an object of class 'arena'!")
}
U = arena$U
V = arena$V
lat = arena$lat
lon = arena$lon
S = arena$S

# Get domain dimension
tmp = dim(S)
m = tmp[1]
n = tmp[2]

num_poly = max(S) +1

# Set the number of subdivisions per time step in such a way that the 
# average jump size is 1/10 of a cell, unless it is specified.
subdiv = max(1, round(t_step * 10 * mean(c(abs(U),abs(V)))))

# Create biomass distribution matrix
N <- matrix(0,nrow=num_poly,ncol=end_t_counter+1)

## Particles are inserted in the C program, according to the inflow rates
## at the boundaries of the domain. This seems the best approximation
## to the introduction of a real biomass into the domain.
## OR: Their prior distribution is given and they are distributed 
## randomly in their polygons according to start_setup.

particles = matrix(0,ncol=2,nrow=num_particles)
results = matrix(0,ncol=12,nrow=num_particles)

# Move particles
if(is.null(start_setup)){
	res_part = .C("particle_tracking_compare_inflow",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),N = as.double(N),as.integer(infl_time),as.integer(infl_poly),as.integer(subdiv),as.double(diffusion))
}
else{
	# start_setup is a vector of biomass distribution over all num_poly polygons
	if(length(start_setup) != num_poly){
		stop("start_setup has wrong length!")
	}
	else{
		# Distribute particles according to start_setup
		# Norm start_setup to 100
		if((sum(start_setup) == 0) | (any(start_setup < 0))){
			stop("The start setup contains no or negative biomass!")
		}
		else{
			start_setup = 100 * start_setup / sum(start_setup)
            }
		res_part = .C("particle_tracking_compare_start_setup",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),N = as.double(N),as.double(start_setup),as.integer(subdiv),as.double(diffusion))
	}
}


# The results contain the following information:
# - The particle positions in part_x and part_y
# - Their positions at each time step in N
particles[,1] = res_part$part_x
particles[,2] = res_part$part_y
N = matrix(res_part$N, nrow=num_poly)

# Norm the number of particle in each time step to 100 (except for the first 
# few, when not all biomass is inserted). They can be different due to 
# stranded particles!
N = N[,1:end_t_counter]
if(is.null(start_setup)){
	if(end_t_counter > infl_time){
		N = N / matrix(pmax(apply(N,2,sum),1),ncol=end_t_counter,nrow=num_poly,byrow=T) * matrix(c(seq(0,100,length=infl_time+1),seq(100,100,length=end_t_counter - infl_time -1)),ncol=end_t_counter,nrow=num_poly,byrow=T)
	}
	else{
		N = N / matrix(pmax(apply(N,2,sum),1),ncol=end_t_counter,nrow=num_poly,byrow=T) * matrix(seq(0,(end_t_counter-1)*100/infl_time,length=end_t_counter),ncol=end_t_counter, nrow=num_poly,byrow=T)
	}
}
else{
	N = N / matrix(pmax(apply(N,2,sum),1),ncol=end_t_counter,nrow=num_poly,byrow=T) * matrix(100, ncol=end_t_counter, nrow=num_poly,byrow=T)
}

if(silent == FALSE){
	# Display N
	BM = S
	dev.new()
	S_index = ifelse(S==0,1,S) # The added index 1 should be filtered out below.
	for (t in 1:end_t_counter){
		BM = ifelse(S==0,0,N[S_index,t])
		BM = matrix(BM,ncol=n,nrow=m)
		readline()
		print(t)
		print(N[,t])
		filled.contour(x=1:n,y=1:m,t(BM),col=terrain.colors(12))
	}
}

# Show the remaining particles and report their number
if(graphics == TRUE){
	dlat = lat[2] - lat[1]
	dlon = lon[2] - lon[1]
	points(particles[,1]*dlon + (lon[1]-dlon),particles[,2]*dlat + (lat[1]-dlat),col="red")
}
if(sum(ifelse(particles[,1] == -10, 0, 1)) > 0){
	warning("Particles left in domain after the time for particle tracking elapsed: ",call. = F)
	warning(sum(ifelse(particles[,1] == -10, 0, 1)),call.=FALSE)
}

result = list(bmdist=N,times=seq(0,(end_t_counter-1)*t_step,by=t_step))
class(result) = "biomass_distribution"
return(result)
}
