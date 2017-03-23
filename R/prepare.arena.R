# Read the flow vectors in east-west (U) and north-south (V) direction 
# from text files, together with latitude and longitude. 
# Additionally, S may be given or is created here.



#' Prepares the arena for an estimation of movement
#' 
#' This function prepares the arena, i.e. the flow field and the polygon
#' structure, for use by the other functions in package BiomassTracking.
#' 
#' The data descriptors \code{file\_U}, \code{file\_V} and \code{file\_S} can
#' come in a number of different formats. If they are filenames, the respective
#' data is loaded from these files. In case the data is already loaded, they
#' can also be data.frames. If they are NULL, a default construction is
#' invoked, where the values \code{n} and \code{m} are needed and the flow is
#' straight in west-east direction.  They can also be matrices giving the
#' respective data directly.
#' 
#' In the usual case, the data is loaded from files and then prepared to be
#' used. The matrices \code{U} and \code{V} are interpolated in every NA
#' position that is not land, while the velocities on land are set to 0. The
#' component perpendicular to land of the sea velocities directly at the
#' land-sea boundary is set to 0 if it points towards land, so that particles
#' wouldn't cross the boundary if an ideal particle tracking algorithm was
#' used.
#' 
#' The data is then rescaled, as the files are supposed to be in m/s format
#' while the particle tracking algorithms act on a longitude/latitude grid.
#' Therefore, the values of 'U' and 'V' are rescaled to degree lon/day and
#' degree lat/day, respectively. Note that this rescaling includes a relative
#' increase of the southern U-velocities compared to the northern U-velocities,
#' because the longitudinal distances are smaller in the South. Also note that
#' the whole algorithm is meant to work in the southern hemisphere and has to
#' be adapted to the Northern Hemisphere if that is wanted.
#' 
#' @param file_U Describes the U data, see details.
#' @param file_V Describes the V data, see details.
#' @param file_S Describes the S data, see details.
#' @param n Zonal dimension of the arena, if needed.
#' @param m Meridional dimension of the arena, if needed.
#' @param na.string The character string denoting NA in the input files.
#' @param eddy If TRUE, an eddy is added in the middle of the arena.
#' @param smooth If TRUE, the flow field is smoothed with a 3x3 moving average.
#' @param noise_sd The standard deviation of Gaussian noise added to the flow
#' field.
#' @return An object of class \code{arena}, describing the arena in which the
#' function is too be used. It is a list containing elements \code{lat},
#' \code{lon}, \code{U}, \code{V} and \code{S} (in that order).  \item{lon,
#' lat}{ Vectors storing the latitude and longitude values of the grid points
#' used.} \item{U, V}{Matrices with the corresponding flow velocities in
#' west-east and south-north direction.} \item{S}{A matrix in which each grid
#' points has an integer number, either giving the polygon it belongs to (if >
#' 0) or stating that this grid point lies on land (if == 0).}
#' @section Warning: This function is written only for the Southern Hemisphere,
#' but can be adapted easily to the Northern Hemisphere.
#' @author Thorsten Lenser
#' @seealso \code{\link{particle.tracking}}, \code{\link{biomass.tracking}},
#' \code{\link{particle.tracking.compare}}, \code{\link{biomass.compare}},
#' \code{\link{plot.arena}}
#' @references TODO: my report
#' @keywords misc
#' @examples
#' 
#' data(Udata)
#' data(Vdata)
#' data(Sdata)
#' arena = prepare.arena(Udata,Vdata,Sdata)
#' plot(arena)
#' 
prepare.arena <- function(file_U=NULL,file_V=NULL,file_S=NULL,n=dim(file_U)[2],m=dim(file_U)[1], na.string="-9999", eddy=FALSE, smooth=FALSE, noise_sd=0){
	# Read the data from the files

	if(is.null(file_U) || is.null(file_V)){
		if(is.null(n) || is.null(m)){
			stop("Specify either U & V or n & m!")
		}
		else{
			# Generate flow field
			#U <- randomflowfield(n,m,0.3)
			#V <- randomflowfield(n,m,0.3)
			#U <- matrix(c(seq(-1,-1,length=n*m/2),seq(1,1,length=n*m/2)),ncol=n)
			#V <- matrix(c(seq(1,1,length=n*m/2),seq(-1,-1,length=n*m/2)),ncol=n, byrow=TRUE)
			U <- matrix(seq(1,1,length=n*m),ncol=n)
			V <- matrix(seq(0,0,length=n*m),ncol=n)
	
			lat = 1:m
			lon = 1:n
		}
	}
	else if(is.matrix(file_U) && is.matrix(file_V)){
		U = file_U
		V = file_V
		m = dim(U)[1]
		n = dim(U)[2]
		if(any(dim(U) != dim(V))){
			stop("If U and V are given as matrices, their dimensions must agree!")
		}
		lat = 1:m
		lon = 1:n
	}
	else{
		if(is.data.frame(file_U)) U = file_U
		else U = read.table(file_U, na.strings = na.string)
		if(is.data.frame(file_V)) V = file_V
		else V = read.table(file_V, na.strings = na.string)
	
		# Extract longitude and latitude information and check 
		# whether they agree for both data sets
		lon = U[,1]
		lat = U[,2]
		if(any(lon != V[,1]) || any(lat != V[,2])){
			stop("Longitude and latitude specification in the two files are not equal!")
		}

		# Keep just the velocity values
		U = U[,3]
		V = V[,3]

		# Get rid of duplicates in lon and lat:
		# set operation union removes duplicates
		lon = unique(lon)
		lat = unique(lat)

		# Determine number of grid points
		n = length(lon)
		m = length(lat)

		# Covert data into matrices
		# U[1,1] should have the smallest lon and lat values
		U = matrix(U,nrow=m,ncol=n,byrow=T)
		V = matrix(V,nrow=m,ncol=n,byrow=T)

		## Rescale the velocity field so that the speed corresponds to the 	
		## distance between station points, i.e. grid points. Since this
		## distance is smaller in the south than in the north, the relative
		## velocities in the south are higher than in the north when equal 
		## in absolute terms.
		## The result is in degrees longitude (for U) or latitude (for V) 
		## per day.
	
		# Convert angles from degrees to radians for rescaling
		lat_rad = lat * pi / 180

		# Turn lat around to fit U
		lat_rad = lat_rad[length(lat):1]
		
		# Scale from s^(-1) to day^(-1) and from km to m, * 86400 / 1000
		U = 86.4 * U
		V = 86.4 * V

		# Scale flow matrices to get degree/day values
		U = U / (111.324 * matrix(cos(lat_rad),nrow=dim(U)[1],ncol=dim(U)[2]))
		V = V / 111.324

		# Scale speed from lon/day or lat/day to grid-point-distance/day
		dlon = lon[2] - lon[1]
		dlat = lat[2] - lat[1]
		U = U / dlon
		V = V / dlat
	}

	# If S is given, 0 denotes land and the polygons are marked
	# 1..max(S). If it is not given, assume that the whole water area
	# is one polygon and the NA values come from land.
	
	if(is.null(file_S)){
		S = array(1,dim=dim(U))
		S = ifelse(is.na(U) | is.na(V),0,S)
	}
	else if(is.data.frame(file_S)){
		S = as.matrix(file_S)
	}
	else if(is.matrix(file_S)){
		S = file_S
		# Determine the extend of land
		S = ifelse(is.na(U) | is.na(V),0,S)
	}
	else if(file_S == "squares"){
		S <- matrix(1,nrow=m,ncol=n)
		S[1:(m/2),1:(n/2)] <- 1
		S[1:(m/2),(n/2 +1):n] <- 2
		S[(m/2 +1):m,1:(n/2)] <- 3
		S[(m/2 +1):m,(n/2 +1):n] <-4
	}
	else if (file_S == "rects"){
		S <- matrix(1,nrow=m,ncol=n)
		S[1:floor(m/3),1:floor(n/2)] <- 1
		S[1:floor(2*m/3),floor(n/2 +1):n] <- 2
		S[floor(m/3 +1):m,1:floor(n/2)] <- 3
		S[floor(2*m/3 +1):m,floor(n/2 +1):n] <- 4
	}
	else if (file_S == "linear"){
		S <- matrix(1,nrow=m,ncol=n)
		S[1:m,1:floor(n/5)] <- 1
		S[1:m,(floor(n/5)+1):floor(2*n/5)] <- 2
		S[1:m,(floor(2*n/5)+1):floor(3*n/5)] <- 3
		S[1:m,(floor(3*n/5)+1):floor(4*n/5)] <- 4
		S[1:m,(floor(4*n/5)+1):n] <- 5
	}
	else{
		S = read.table(file_S, na.strings = na.string)
		# convert to pure matrix, there must be a better way...
		S = array(as.integer(as.matrix(S)),dim=dim(S)) 
	}

	## Now interpolates the U and V data in regions where the 	
	## topography is too shallow to get measured results.

	# Interpolate any NA values which are in shallow water, not on land.
	while(any(is.na(U) & (S != 0))){
		for(i in 1:m){
			for(j in 1:n){
				if(is.na(U[i,j]) && (S[i,j] != 0)){
					tmp = 0
					count = 0
					if((i > 1) && (! is.na(U[i-1,j]))){
						tmp = tmp + U[i-1,j]
						count = count + 1
					}
					if((i < m) && (! is.na(U[i+1,j]))){
						tmp = tmp + U[i+1,j]
						count = count + 1
					}
					if((j > 1) && (! is.na(U[i,j-1]))){
						tmp = tmp + U[i,j-1]
						count = count + 1
					}
					if((j < n) && (! is.na(U[i,j+1]))){
						tmp = tmp + U[i,j+1]
						count = count + 1
					}
					if(count > 0) U[i,j] = tmp/count
				}
			}
		}
	}
	while(any(is.na(V) & (S != 0))){
		for(i in 1:m){
			for(j in 1:n){
				if(is.na(V[i,j]) && (S[i,j] != 0)){
					tmp = 0
					count = 0
					if((i > 1) && (! is.na(V[i-1,j]))){
						tmp = tmp + V[i-1,j]
						count = count + 1
					}
					if((i < m) && (! is.na(V[i+1,j]))){
						tmp = tmp + V[i+1,j]
						count = count + 1
					}
					if((j > 1) && (! is.na(V[i,j-1]))){
						tmp = tmp + V[i,j-1]
						count = count + 1
					}
					if((j < n) && (! is.na(V[i,j+1]))){
						tmp = tmp + V[i,j+1]
						count = count + 1
					}
					if(count > 0) V[i,j] = tmp/count
				}
			}
		}
	}

	## This part modifies the matrices U and V in such a way that there is 
	## a no-pass boundary condition in place at the land-sea edge, i.e.
	## the perpendicular flow component from sea onto land is set to 0.
	## Note: Unless another particle tracking algorithm is used, particles
	## will still be able to cross the land-sea boundary, due to the fixed time 
	## step(?). In this case, they will be treated as if they left the domain.
	## Land cells are marked with the maximum entry in S.

	for (i in 1:m){
		for (j in 1:n){
			if(S[i,j] == 0){
				if((i > 1) && (S[i-1,j] != 0)){
					if(V[i-1,j] > 0) V[i-1,j] = 0
				}
				if((i < m) && (S[i+1,j] != 0)){
					if(V[i+1,j] < 0) V[i+1,j] = 0
				}
				if((j > 1) && (S[i,j-1] != 0)){
					if(U[i,j-1] > 0) U[i,j-1] = 0
				}
				if((j < n) && (S[i,j+1] != 0)){
					if(U[i,j+1] < 0) U[i,j+1] = 0
				}

				# Set all the flow on land to 0
				U[i,j] = 0
				V[i,j] = 0
			}
		}
	}				

	# Shallow waters and the coastal regions are filled, fill the rest with
	# zeros.
	U = ifelse(is.na(U),0,U)
	V = ifelse(is.na(V),0,V)

	# Put an eddy into middle field
	if(eddy == TRUE){
		for (i in 1:(m/3)){
			for (j in 1:(n/3)){
				U[i+(m/3),j+(n/3)] <- -0.5*(i - m/6)
				V[i+(m/3),j+(n/3)] <- 0.5*(j - n/6)
			}
		}
	}

	# Smooth the flow field
	if(smooth == TRUE){
		for (i in 2:(m-1)){
			for (j in 2:(n-1)){
				U[i,j] <- (U[i-1,j-1]+U[i-1,j]+U[i-1,j+1]+U[i,j-1]+U[i,j]+U[i,j+1]+U[i+1,j-1]+U[i+1,j]+U[i+1,j+1])/9
				V[i,j] <- (V[i-1,j-1]+V[i-1,j]+V[i-1,j+1]+V[i,j-1]+V[i,j]+V[i,j+1]+V[i+1,j-1]+V[i+1,j]+V[i+1,j+1])/9
			}
		}
	}

	# Add a bit of noise
	if(noise_sd > 0){
		U = U + rnorm(n*m,0,noise_sd)
		V = V + rnorm(n*m,0,noise_sd)
	}


	# Export the values
	result = list(lon=lon,lat=lat,U=U,V=V,S=S)
	class(result) = "arena"
	return(result)
}
