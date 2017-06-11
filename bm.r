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
are <- prepare.arena(Udata, Vdata, Sdata)
pp <- pm(are, 40, 5, diffusion = 0.06, graphics = TRUE, from_boundary = FALSE, log = FALSE)
pm <- function(arena,num_particles,no_move_window=1000,end_t_counter=NULL,t_step=NULL,diffusion=0,no_age_classes="FD",min_nac=100,max_nac=500,delta=NULL,graphics=FALSE,log=FALSE,from_boundary=TRUE){
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
  num_poly = max(S)+1
  
  # Get domain dimension
  m = dim(S)[1]
  n = dim(S)[2]
  if(from_boundary == TRUE){
    # Particles are entered in the C program
    particles = matrix(0,ncol=2,nrow=num_particles)
  }
  else{
    
    # Insert num_particles particles randomly all over the domain
    particles = matrix(c(runif(num_particles,min=0.501,max=(n + 0.499)),runif(num_particles,min=0.501,max=(m + 0.499))),ncol=2)
    print(str(particles))
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
  
res
}
