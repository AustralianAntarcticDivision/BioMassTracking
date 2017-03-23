## particle.tracking.R
## author: Thorsten Lenser (thorsten.lenser@gmx.net)

## The particle tracking algorithm with inflow and outflow at the 
## domain boundaries. Estimates the movement tensors after tracking
## the paths of a given number of particles.
## Uses C functions to do the hard work.

## There are as many moves as there are polygons (including the outside).
## The time, origin and direction of of all moves are recorded and the 
## estimates are made from these.

## Land (or generally regions with NA data points) is treated as outside! 
## It is marked with 0 in matrix S, but is then treated as polyon max(S)+1 !!!

## For help, see the documentation.



#' Estimate a 'movement kernel' describing the movement of biomass.
#' 
#' A \code{movement kernel} is meant to summarise the movement of biomass
#' across polygon-shaped areas in the ocean. It is estimated using a 4-th order
#' Runge-Kutta particle tracking algorithm, which records the moves of passive
#' numerical drifters from polygon to polygon.
#' 
#' The two main bits of work in this function, the Runge-Kutta simulation and
#' the estimation of the \code{movement\_ kernel}, have been implemented in the
#' C dlls \code{particle\_tracking} and \code{estimate\_P}. Even so, the
#' runtime of this function can be quite high when using large numbers of
#' particles and/or a large domain. It is important to note that in the present
#' form, particles can get caught in 'voids' in the flow field, i.e. in places
#' where their velocity is 0. If this happens for too many particles the
#' estimation should be repeated with a higher values of diffusion.  A warnings
#' message is given at the end of the simulation, indicating the number of
#' particles left in the domain.
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
#' The default values for the time step in the Runge-Kutta scheme is choosen in
#' such a way that a particle will on average move for 1/10 of the distance
#' between two grid points in each time step. Although a variable time step
#' would be preferable, this yields an accurate estimation in most cases.
#' 
#' Diffusion is included by adding the product of a normal random number and
#' the mean velocity in the either direction to the calculated velocity in that
#' direction. The diffusion coefficient gives the standard variation of the
#' normal random number. Diffusion is needed to get particles out of 'voids',
#' i.e. places in the flow field where the adjacent velocities cancel out and
#' the particles get stuck. It is also physically justifiable as the flow field
#' itself characterises advective motion only.
#' 
#' The main loss of precision in the replication and prediction of movement
#' stems from the management of time, which is realised as a system of 'age
#' classes' in the movement kernel. A particle entering a polygon is in age
#' class 1 and then moves up through the age classes. Each age class has
#' different probabilities of movement attached to it, so that the movement of
#' a biomass depends on the time since it entered its current polygon. Since
#' the choice of the number of age classes and their width \code{delta} can
#' have a profound effect, the user can choose whether to trust one of the
#' built-in methods (\code{Sturges}, \code{FD}, \code{scott}) which are usually
#' used to choose the number of classes in a histogram, or specificy his or her
#' own number of age classes. The default uses the Freedman-Diaconis rule
#' (\code{FD}), which usually produces the highest number of age-classes and
#' therefore the best precision. Unless the performance of the function
#' \code{biomass\_tracking} is totally unsatisfactory due to the high number of
#' age-classes, less age-classes should not be used, as the loss of precision
#' can be quite severe.
#' 
#' @param arena An object of class \code{arena}, see details.
#' @param num_particles The number of particles used for the estimation.
#' @param no_move_window Number of time steps in which no move between polygons
#' has to occur in order for the particle tracking to stop. Use 0 to disable.
#' @param end_t_counter The maximal number of time steps used in the
#' Runge-Kutta scheme. Has to be specified if log==TRUE.
#' @param t_step The time step used in the Runge-Kutta scheme. See details for
#' the default value.
#' @param diffusion A diffusion coefficient, see details.
#' @param no_age_classes The number of age classes used in the \code{movement
#' kernel}, takes \code{Sturges}, \code{FD}, \code{scott} or an integer number.
#' See details.
#' @param min_nac The minimum number of age classes.
#' @param max_nac The maximum number of age classes.
#' @param delta The width of an age class used in the movement kernel, chosen
#' automatically if NULL.
#' @param graphics If TRUE, the remaining particles at the end of the
#' Runge-Kutta simulation are shown.
#' @param log If TRUE, the path of the particles is saved to a global variable
#' 'log' and displayed if graphics==TRUE. Only use for a small number of
#' particles!
#' @param from_boundary If TRUE, particles are inserted from the boundary,
#' otherwise they are uniformly distributed over the domain. This might be
#' advantageous for plotting purposes if log=TRUE.
#' @return An object of class \code{movement\_kernel}, which is a list of the
#' following items: \item{P}{A list of the movement tensors for each polygon,
#' where \code{P[[poly]][i,j,k]} describes the probability of a particle
#' leaving polygon \code{poly} towards polygon \code{j} in the case that it
#' entered from polygon \code{i} and is in age class \code{k}.} \item{nk}{A
#' list containing the neighbours for each polygon. This is used to translate
#' between the local neighbourhood numbers (used in \code{P}) to the global
#' numbers of the polygons.} \item{no_age_classes}{The number of age classes
#' used in the movement kernel.} \item{delta}{Width of the age classes used in
#' the movement kernel.}
#' @author Thorsten Lenser
#' @seealso \code{\link{biomass.tracking}},
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
#' @export  
particle.tracking <- function(arena,num_particles,no_move_window=1000,end_t_counter=NULL,t_step=NULL,diffusion=0,no_age_classes="FD",min_nac=100,max_nac=500,delta=NULL,graphics=FALSE,log=FALSE,from_boundary=TRUE){
  times1 <- times_save <- origin_save <- cur_poly_save <- dest_save <- NULL
  if(is.null(t_step)){
    t_step = 1/(10*mean(c(abs(arena$U),abs(arena$V))))
  }
  
  if(class(arena) != "arena"){
    stop("First argument has to be an object of class 'arena'!")
  }
  U = arena$U
  V = arena$V
  lat = arena$lat
  lon = arena$lon
  S = arena$S
  
  if(is.null(end_t_counter) && (log == TRUE)){
    stop("If log==TRUE, end_t_counter has to be given!")
  }
  
  # Get domain dimension
  m = dim(S)[1]
  n = dim(S)[2]
  
  # Get number of polygons, including outside and land, which form one polygon
  num_poly = max(S)+1
  
  # Get local neighbourhood numbers for each polygon
  nk = list()
  for(poly in 1:num_poly){
    nk = c(nk,list(numeric(0)))
  }
  for(i in 1:m){
    for(j in 1:n){
      if(S[i,j] != 0){
        if(i < m){
          if(S[i+1,j] != 0) nk[[S[i,j]]] = union(nk[[S[i,j]]],S[i+1,j])
          else nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)	
        }
        else{
          nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        if(i > 1){
          if(S[i-1,j] != 0) nk[[S[i,j]]] = union(nk[[S[i,j]]],S[i-1,j])		
          else nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        else{
          nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        if(j < n){
          if(S[i,j+1] != 0) nk[[S[i,j]]] = union(nk[[S[i,j]]],S[i,j+1])		
          else nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        else{
          nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        if(j > 1){
          if(S[i,j-1] != 0) nk[[S[i,j]]] = union(nk[[S[i,j]]],S[i,j-1])		
          else nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
        else{
          nk[[S[i,j]]] = union(nk[[S[i,j]]],num_poly)
        }
      }
      else{
        if(i < m){
          if(S[i+1,j] != 0) nk[[num_poly]] = union(nk[[num_poly]],S[i+1,j])	
        }
        else{
          nk[[num_poly]] = union(nk[[num_poly]],num_poly)
        }
        if(i > 1){
          if(S[i-1,j] != 0) nk[[num_poly]] = union(nk[[num_poly]],S[i-1,j])		
        }
        else{
          nk[[num_poly]] = union(nk[[num_poly]],num_poly)
        }
        if(j < n){
          if(S[i,j+1] != 0) nk[[num_poly]] = union(nk[[num_poly]],S[i,j+1])		
        }
        else{
          nk[[num_poly]] = union(nk[[num_poly]],num_poly)
        }
        if(j > 1){
          if(S[i,j-1] != 0) nk[[num_poly]] = union(nk[[num_poly]],S[i,j-1])		
        }
        else{
          nk[[poly]] = union(nk[[num_poly]],num_poly)
        }
      }	
    }
  }
  for(i in 1:m){
    nk[[num_poly]] = union(nk[[num_poly]],S[i,1])
    nk[[num_poly]] = union(nk[[num_poly]],S[i,n])
  }
  for(j in 1:n){
    nk[[num_poly]] = union(nk[[num_poly]],S[1,j])
    nk[[num_poly]] = union(nk[[num_poly]],S[m,j])
  }
  
  # Add the polygons as their own neighbour to accommodate native biomass
  for(poly in 1:num_poly){
    nk[[poly]] = union(nk[[poly]],poly)
  }
  nk[[num_poly]] = setdiff(nk[[num_poly]],0)
  
  
  if(from_boundary == TRUE){
    # Particles are entered in the C program
    particles = matrix(0,ncol=2,nrow=num_particles)
  }
  else{
    # Insert num_particles particles randomly all over the domain
    particles = matrix(c(runif(num_particles,min=0.501,max=(n + 0.499)),runif(num_particles,min=0.501,max=(m + 0.499))),ncol=2)
    
    # Regular distribution of particles, see handwriting
    #part_cols = floor(sqrt(num_particles*n/m))
    #part_rows = floor(m*part_cols/n)
    #dist_x = (n+1) / part_cols
    #dist_y = (m+1) / part_rows
    #x_part = seq(0.5 + (dist_x/2),((n+0.5) - (dist_x/2)),length=part_cols)
    #y_part = seq(0.5 + (dist_y/2),((m+0.5) - (dist_y/2)),length=part_rows)
    #num_particles = part_cols * part_rows
    #particles = matrix(0,ncol=2,nrow=num_particles)
    
    #for (i in 1:part_rows){
    # 	for (j in 1:part_cols){
    #		particles[(i-1)*part_cols +j,1] = x_part[j]
    #		particles[(i-1)*part_cols +j,2] = y_part[i]
    #	}
    #}
  }
  
  # See if a maximum number of iterations has been specified.
  if(is.null(end_t_counter)){
    end_t_counter = 0
  }
  
  # Create log data storage if needed
  if(log == TRUE){
    logdata = matrix(0,ncol=num_particles * end_t_counter, nrow=2)
  }
  
  # Move particles & record the movements
  if(log == FALSE){
    if(from_boundary == TRUE){
      res = .C("particle_tracking",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(no_move_window),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),res_t = as.double(numeric(num_particles*num_poly)),res_S1 = as.integer(integer(num_particles*num_poly)),res_S2 = as.integer(integer(num_particles*num_poly)),res_S3 = as.integer(numeric(num_particles*num_poly)),as.double(diffusion),as.integer(0),as.double(0),as.integer(1))
    }
    else{
      res = .C("particle_tracking",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(no_move_window),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),res_t = as.double(numeric(num_particles*num_poly)),res_S1 = as.integer(integer(num_particles*num_poly)),res_S2 = as.integer(integer(num_particles*num_poly)),res_S3 = as.integer(numeric(num_particles*num_poly)),as.double(diffusion),as.integer(0),as.double(0),as.integer(0))
    }
  }
  else{
    no_move_window = 0
    if(from_boundary == TRUE){
      res = .C("particle_tracking",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(no_move_window),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),res_t = as.double(numeric(num_particles*num_poly)),res_S1 = as.integer(integer(num_particles*num_poly)),res_S2 = as.integer(integer(num_particles*num_poly)),res_S3 = as.integer(numeric(num_particles*num_poly)),as.double(diffusion),as.integer(1),logdata=as.double(logdata),as.integer(1))
    }
    else{
      res = .C("particle_tracking",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(no_move_window),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),res_t = as.double(numeric(num_particles*num_poly)),res_S1 = as.integer(integer(num_particles*num_poly)),res_S2 = as.integer(integer(num_particles*num_poly)),res_S3 = as.integer(numeric(num_particles*num_poly)),as.double(diffusion),as.integer(1),logdata=as.double(logdata),as.integer(0))
    }
  }
  
  # The results contain the following information:
  # - The final particle positions in part_x and part_y
  # - Time of the movements in res_t
  # - Where the moving particles came from in res_S1
  # - The current polygon of the moving particles in res_S2
  # - Destination of the movement in res_S3
  #
  # Each of these if has all moves of one 'order' in one block,
  # i.e. the first move is in the first block of 'num_particles' entries.
  #
  # Store this information in particles, times, origin, cur_poly and dest
  particles[,1] = res$part_x
  particles[,2] = res$part_y
  
  # Calculate the relative movement times from the absolute ones
  for(t in num_poly:2){
    res$res_t[((t-1)*num_particles+1):(t*num_particles)] = ifelse(res$res_t[((t-1)*num_particles+1):(t*num_particles)] > 0, ifelse(res$res_S3[((t-1)*num_particles+1):(t*num_particles)] != 0, res$res_t[((t-1)*num_particles+1):(t*num_particles)] - res$res_t[((t-2)*num_particles+1):((t-1)*num_particles)], res$res_t[((t-1)*num_particles+1):(t*num_particles)]), 0)
  }
  
  # Store the results in separate variables
  times = res$res_t
  origin = res$res_S1
  cur_poly = res$res_S2
  dest = res$res_S3
  
  # Get logdata if used
  if(log==TRUE){
    logdata = res$logdata
  }
  
  rm(res)
  
  # Throw out all the records where no move occured (because the
  # particle got stuck before).
  # To make sure we keep the 1st record of a particles when it got stuck,
  # use times==0, and not dest==0. Otherwise the estimation does not work properly.
  if(any(dest == 0)) some_did_not_move = TRUE
  else some_did_not_move = FALSE
  origin = as.integer(na.omit(ifelse(times==0,NA,origin)))
  cur_poly = as.integer(na.omit(ifelse(times==0,NA,cur_poly)))
  dest = as.integer(na.omit(ifelse(times==0,NA,dest)))
  times = as.numeric(na.omit(ifelse(times==0,NA,times)))
  
  # Numerical errors could have been introduced to times, due to 
  # substraction. Keep this in mind when debugging!
  
  # Distance between two grid points, for graphics
  dlat = lat[2] - lat[1]
  dlon = lon[2] - lon[1]
  
  # Show the particles that are left in the domain, print how many are left
  if(graphics == TRUE){
    points(particles[,1]*dlon + (lon[1]-dlon),particles[,2]*dlat + (lat[1]-dlat),col="red")
  }
  if(sum(ifelse(particles[,1] == -10, 0, 1)) > 0){
    warning("Particles left in domain after the time for parameter estimation elapsed: ",call. = F)
    warning(sum(ifelse(particles[,1] == -10, 0, 1)),call.=FALSE)
  }
  
  # Now estimate the movement of native biomass that "grows" inside a polygon,
  # i.e. did not enter from anywhere else. 
  particles = matrix(c(runif(num_particles,min=0.501,max=(n + 0.499)),runif(num_particles,min=0.501,max=(m + 0.499))),ncol=2)
  res = .C("particle_tracking_native_biomass",as.double(U),as.double(V),as.integer(n),as.integer(m),as.integer(num_particles),as.integer(S),as.double(t_step),as.integer(no_move_window),as.integer(end_t_counter),part_x = as.double(particles[,1]),part_y = as.double(particles[,2]),res_t = as.double(numeric(num_particles)),res_S2 = as.integer(integer(num_particles)),res_S3 = as.integer(numeric(num_particles)),as.double(diffusion))
  monb_times = res$res_t
  monb_cur_poly = res$res_S2
  monb_dest = res$res_S3
  
  rm(res)
  
  # Throw out all the records where no move occured (because the
  # particle got stuck before)
  if(any(monb_dest == 0)) some_did_not_move = TRUE
  monb_cur_poly = as.integer(na.omit(ifelse(monb_times==0,NA,monb_cur_poly)))
  monb_dest = as.integer(na.omit(ifelse(monb_times==0,NA,monb_dest)))
  monb_times = as.numeric(na.omit(ifelse(monb_times==0,NA,monb_times)))
  
  # Combine both sets of movement data
  times = c(times, monb_times)
  origin = c(origin, monb_cur_poly)
  cur_poly = c(cur_poly, monb_cur_poly)
  dest = c(dest, monb_dest)
  times1 <<- times
  # 'Cut off' the times data at the point where the quantile difference between two neighbouring quantiles
  # (in 0.01 increment) is greatest. This point is likely to be the point where the movement more or less 
  # stopped. Do this only if some particles got retained.
  if(some_did_not_move == T){
    max_t_index = which.max(quantile(times,seq(0.01,1,by=0.01)) - quantile(times,seq(0,0.99,by=0.01)))
    max_t_index = min(100, max_t_index)
    max_t_index = max(2, max_t_index)
    
    # Use the average of the three points around max_t_index-1 (the jump
    # occurs before max_t_index). This is needed in cases of very 'spikey' data.
    times = pmin(times, sum(quantile(times, c(max_t_index-2, max_t_index-1 , max_t_index)/ 100))/3 )
  }
  times = pmin(times,5000)
  
  # Estimate the number of age classes  and delta
  if(no_age_classes == "Sturges"){
    no_age_classes = nclass.Sturges(times)
  }
  else if(no_age_classes == "scott"){
    no_age_classes = nclass.scott(times)
  }
  else if(no_age_classes == "FD"){
    no_age_classes = nclass.FD(times)
  }
  # Sometimes this does not work because all the movement values are equal,
  # in which case the estimation fails, or are very near together, in which case
  # the default 'FD' yields a very high number of age-classes.
  # Therefore, we use a minimum and a maximum number of age-classes
  if(is.na(no_age_classes) | is.nan(no_age_classes) | (no_age_classes < min_nac)){
    no_age_classes = min_nac
  }
  if(no_age_classes > max_nac) no_age_classes = max_nac
  
  if((some_did_not_move==TRUE) && is.null(delta)){
    # Choose delta in such a way that the last age class contains only the
    # particles which did not move, which have the maximal time stamp.
    delta = max((max(ifelse(times == max(times), NA, times), na.rm=T) + max(times))/2, max(times)/(1 + 1/(2*no_age_classes -2))) / (no_age_classes -1)
  }
  else if(is.null(delta)){
    # Choose delta in such a way that the last age class contains no 
    # no particles, so that nothing is retained.
    delta = (max(times) + 1) / (no_age_classes -1)
  }
  
  
  # Store data for debugging
  times_save <<- times
  origin_save <<- origin
  cur_poly_save <<- cur_poly
  dest_save <<- dest
  
  # Process the log data and display it if wanted
  # Change -10 values to NA, so it is clear that these particles left the domain.
  if(log == TRUE){
    logdata = ifelse(logdata == -10, NA, logdata)
    logdata = matrix(logdata, nrow=2)
    log <<- logdata
    if(graphics == TRUE){
      for(i in 1:num_particles){
        lines(x = (logdata[1,(0:(end_t_counter-1)) * num_particles +i])*dlon + (lon[1]-dlon), y = (logdata[2,(0:(end_t_counter-1)) * num_particles +i])*dlat + (lat[1]-dlat))
      }
    }
  }
  
  
  # Create movement probability tensors, one for each polygon, 
  # with indices origin,dest,t (in that order!)
  P = list()
  for(poly in 1:(num_poly-1)){
    P = c(P,list(array(0, dim = c(length(nk[[poly]]),length(nk[[poly]]),no_age_classes))))
  }
  
  # Estimate movement tensors.
  for (poly in 1:(num_poly-1)){
    P[[poly]] = .C("estimate_P", P=as.double(P[[poly]]), as.integer(nk[[poly]]), as.integer(na.omit(ifelse(cur_poly==poly,origin,NA))), as.integer(na.omit(ifelse(cur_poly==poly,dest,NA))), as.double(na.omit(ifelse(cur_poly==poly,times,NA))), as.integer(no_age_classes), as.integer(length(as.integer(na.omit(ifelse(cur_poly==poly,origin,NA))))), as.integer(length(nk[[poly]])), as.double(delta), as.integer(1))$P
    P[[poly]] = array(P[[poly]],dim = c(length(nk[[poly]]), length(nk[[poly]]),no_age_classes))
  }
  
  if(some_did_not_move == TRUE){
    result = list(P=P,nk=nk,no_age_classes=no_age_classes,delta=delta,retention=TRUE,diffusion=diffusion)
  }
  else{
    result = list(P=P,nk=nk,no_age_classes=no_age_classes,delta=delta,retention=FALSE,diffusion=diffusion)
    
  }
  class(result) = "movement_kernel"
  return(result)
}

