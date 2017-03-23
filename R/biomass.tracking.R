## Track the movement of biomass through the area, using the movement tensors
## estimated by "particle_tracking".

## Note that this model ticks in time steps of delta, i.e. the width
## of the age classes. The required time points will then be linearly
## interpolated between the calculated ones.

## If silent==TRUE, no text or graphical output will apear.
## S is needed of silent == FALSE!



#' Calculate the movement of biomass from movement tensors
#' 
#' This function calculates the movement of biomass between given polygons as
#' determined by the movement parameters. It uses the movement tensor
#' calculated by a numerical tracer experiment implemented by the function
#' \code{particle.tracking}.
#' 
#' Two different starting conditions for the flow of biomass can be given: A
#' starting distribution of biomass can be specified or it can flow into the
#' area through the specified inflow polygon over a given number of time steps.
#' In this case, time steps refers to the differences between the first two
#' entries in the \code{times} vector.
#' 
#' An object of class \code{movement\_kernel} is a list of the following items:
#' 
#' P: A list of the movement tensors for each polygon, where P[[poly]][i,j,k]
#' describes the probability of a particle leaving polygon \code{poly} towards
#' polygon \code{j} in the case that it entered from polygon \code{i} and is in
#' age class \code{k}. nk: A list containing the neighbours for each polygon.
#' This is used to translate between the local neighbourhood numbers (used in
#' P) to the global numbers of the polygons. no\_age\_classes: The number of
#' age classes used in the movement kernel. delta: Width of the age classes
#' used in the movement kernel. retention: Logical to indicate whether
#' retention in the final age class has been enabled or not.
#' 
#' @param mk An object of class \code{movement kernel}, see details.
#' @param times A vector giving the time-points at which the output is wanted.
#' @param start_setup A vector giving a starting biomass distribution. Specify
#' either this or the next two arguments.
#' @param infl_poly The biomass flows in from this polygon. Has to have
#' positive inflow from the boundary.
#' @param infl_time The number of time steps over which the biomass is to be
#' inserted. Only used if start\_setup==NULL.
#' @param insert_t_step The width of the time steps over which the biomass is
#' to be inserted. Only used if start\_setup==NULL.
#' @param silent If FALSE, the movement is shown in text and graphic and the
#' user has to hit 'enter' after every time step.
#' @param S Matrix describing the polygon structure, needed if silent == FALSE.
#' @return A matrix of dimension 'number of polygons' x 'number of time steps
#' requested'. Each column contains the biomass distribution in the respective
#' time step.
#' @author Thorsten Lenser
#' @seealso \code{\link{particle.tracking}},
#' \code{\link{particle.tracking.compare}}, \code{\link{biomass.compare}},
#' \code{\link{prepare.arena}}
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
biomass.tracking <- function(mk,times,start_setup=NULL,infl_poly=1,infl_time=10,insert_t_step=NULL,silent=TRUE,S=NULL){

if(class(mk) != "movement_kernel"){
	stop("First argument has to be of class 'movement_kernel'!")
}

# Check if the times are stricty increasing
for(i in 2:length(times)){
	if(!(times[i] > times[i-1])) stop("The 'times' have to be strictly increasing!")
}

P = mk$P
nk = mk$nk
delta = mk$delta
no_age_classes = mk$no_age_classes
num_poly = length(P) +1  # including "outside"
if(is.null(S) == FALSE){
	m = dim(S)[1]
	n = dim(S)[2]
	S = ifelse(S==0,max(S)+1,S)
}

# Prepare an array with the P data from the movement kernel
P_array = numeric(0)
for(i in 1:length(P)){
	P_array = c(P_array,P[[i]])
}

# Prepare an array containing the numbers of neighbours 
num_neighbours = integer(0)
for(i in 1:num_poly){
	num_neighbours = c(num_neighbours,length(nk[[i]]))
}

# Prepare an array with nk in it
nk_array = integer(0)
for(i in 1:num_poly){
	nk_array = c(nk_array,nk[[i]]);
}

# Choose a number of time steps for insertion of none is given
if(is.null(insert_t_step)){
	if(length(times) > 1) insert_t_step = times[2] - times[1]
	else insert_t_step = 1
}


# Run the delta-time steps
if(is.null(start_setup)){
	N_delta = .C("biomass_tracking_inflow", N_delta = as.double(matrix(0,nrow=num_poly,ncol=as.integer(max(times) / delta +3))), as.double(max(times)+(delta/2)),as.double(P_array),as.integer(num_neighbours),as.integer(num_poly),as.integer(nk_array),as.integer(no_age_classes),as.double(delta),as.integer(infl_time),as.double(insert_t_step),as.integer(infl_poly))$N_delta
}
else{
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
		N_delta = .C("biomass_tracking_start_setup", N_delta = as.double(matrix(0,nrow=num_poly,ncol=as.integer(max(times) / delta +3))), as.double(max(times)+(delta/2)),as.double(P_array),as.integer(num_neighbours),as.integer(num_poly),as.integer(nk_array),as.integer(no_age_classes),as.double(delta),as.double(start_setup))$N_delta
	}
}
N_delta = matrix(N_delta,nrow=num_poly)

# Interpolate the calculated time steps to the desired time steps
# Assume the calculated time steps took place in the middle of the age-class time-step,
# i.e. at delta/2, 3*delta/2, ...
delta_counter = 1  # count the delta time steps
N <- matrix(0,nrow=num_poly,ncol=length(times))
for(counter in 1:length(times)){
	exact_time = times[counter]
	while((delta_counter * delta - (delta/2)) <= exact_time){
		delta_counter = delta_counter + 1
	}

	if(delta_counter == 1){
		N[,counter] <- N_delta[,delta_counter] + (exact_time  / delta) * (N_delta[,delta_counter+1] - N_delta[,delta_counter])
	}
	else{
		N[,counter] <- N_delta[,delta_counter] + ((exact_time - (delta_counter-1) * delta + delta/2) / delta) * (N_delta[,delta_counter+1] - N_delta[,delta_counter])
	}
	if(silent == FALSE){
		BM = N[as.integer(S),counter]
		BM = matrix(BM,ncol=n,nrow=m)
		#print(counter)
		print(N[,counter])
		image(x=1:n,y=1:m,t(BM),col=terrain.colors(12))
		readline()
	}
}

result = list(bmdist=N,times=times)
class(result) = "biomass_distribution"
return(result)
}
