## This function creates a random mxn matrix meant to resemble one coordinate
## in a flow field.
## It starts by choosing a starting value randomly between -1 and 1 and then 
## performs a random walk on the differences in horizontal and
## vertical direction.

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
