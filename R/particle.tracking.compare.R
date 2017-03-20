## This version of the particle tracking algorithms is not used to estimate
## parameters, but to compare against a coarser model of biomass movement.
##
## If set to silent, no text or graphical output will appear.

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
subdiv = max(1, round(t_step * 10 * mean(abs(U),abs(V))))

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
	windows()
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
