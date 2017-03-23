## This function creates a random mxn matrix meant to resemble one coordinate
## in a flow field.
## It starts by choosing a starting value randomly between -1 and 1 and then 
## performs a random walk on the differences in horizontal and
## vertical direction.



#' Generate one direction of a random flow field
#' 
#' Generate a matrix meant to be used as one direction of a random flow field.
#' 
#' This function proceeds in the following way: A starting value, either 1 or
#' -1, is drawn at random. This is assigned to the element in the top-left
#' corner of the matrix. The other values are then calculated as the sum of a
#' normal random variable with standard deviaton \code{std\_dev} and their left
#' neighbour (top row) or the mean of their left and their upper neighbour
#' (rest).
#' 
#' @param n Zonal (x) extent of the flow field.
#' @param m Meridional (y) extent of the flow field.
#' @param std_dev Standard deviation used for the random walk, see details.
#' @return A matrix meant to be used as one direction of a random flow field.
#' @author Thorsten Lenser
#' @seealso \code{\link{prepare.arena}}
#' @keywords misc
#' @examples
#' 
#' U = randomflowfield(30,30,0.3)
#' V = randomflowfield(30,30,0.3)
#' quiver(U,V)
#' 
randomflowfield <- function(n,m,std_dev){
	# Starting value
	x = runif(1, min=-1,max=1)

	# Ensure a norm of 1 for the starting value
	if(x != 0){
		x = x / abs(x)
	}
	else{
		x = 1
	}
	
	# Get the rest by random walking
	for (j in 2:m){
		x[j] = x[j-1] + rnorm(n=1,sd=0.2)
	}
	for (i in 2:n){
		for (j in 1:m){
			x[(i-1)*m + j] = (x[(i-1)*m + j -1] + x[(i-2)*m + j])/2 + rnorm(n=1,sd = std_dev)
		}
	}

	A = matrix(x, nrow = m, byrow=TRUE)
	return(A)
}
