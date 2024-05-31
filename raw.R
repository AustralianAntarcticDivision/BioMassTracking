library(BiomassTracking)
?particle.tracking



res = .C("particle_tracking", 
         as.double(U),
         as.double(V),
         as.integer(n), ## nrows
         as.integer(m), ## ncols
         as.integer(num_particles),
         
         as.integer(S), ## polygon ID
         as.double(t_step),
         as.integer(no_move_window), ## number of time steps in which no move between polygons has to occur in order for the particle tracking to stop. Use 0 to disable.
         as.integer(end_t_counter),
         part_x = as.double(particles[,1]),
         part_y = as.double(particles[,2]),
         res_t = as.double(numeric(num_particles*num_poly)),
         res_S1 = as.integer(integer(num_particles*num_poly)),
         res_S2 = as.integer(integer(num_particles*num_poly)),
         res_S3 = as.integer(numeric(num_particles*num_poly)),
         as.double(diffusion),
         as.integer(0),as.double(0),as.integer(0))



data(Udata)
data(Vdata)
data(Sdata)

arena = prepare.arena(Udata,Vdata,Sdata)
plot(arena)

# Use more particles if realistic results are needed
#res <- particle.tracking0(arena,4000,5000,diffusion=0.6,graphics=FALSE, from_boundary = T)
#runge_kutta(double *x, double *y, double *U, double *V, int *m, int *n, double *t_step)

num_particles <- 1000
library(raadtools)
x <- readcurr(xylim = extent(120, 180, -73, -20))
x[is.na(x)] <- 0
arena <- list(U = as.matrix(x[[1]])[nrow(x):1, ], 
              V = as.matrix(x[[2]])[nrow(x):1, ])
m <- nrow(arena$U)
n <- ncol(arena$U)
particles0 <-  
  matrix(c(runif(num_particles,
                 min=0.501,max=(n + 0.499)),
           runif(num_particles,min=0.501,max=(m + 0.499))),ncol=2)

t_step <- 1/(10*mean(c(abs(arena$U),abs(arena$V))))
#t_step <- 40
particles <- particles0
plot(particles, col = rgb(0.7, 0.7, 0.7))

library(raadtools)

run_rk <- function(particles)  {
  l <- vector("list", nrow(particles))
  for (i in seq_len(nrow(particles))) {
  rk <- .C("runge_kutta", as.double(particles[i,1]), as.double(particles[i,2]), 
            as.double(arena$U), as.double(arena$V), 
            as.integer(nrow(arena$U)), as.integer(ncol(arena$U)), 
            as.double(t_step))
  l[[i]] <- cbind(rk[[1]], rk[[2]])
  }
  do.call(rbind, l)
}
p <- run_rk(particles)
plot(p, col = "lightgrey")

for (j in 1:400 ) {
  p <- run_rk(p)
  points(p, pch = ".")
}
