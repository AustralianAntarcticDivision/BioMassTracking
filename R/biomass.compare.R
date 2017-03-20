# Compare the two runs of the models

biomass.compare <- function(N1,N2,graphics=TRUE,legend=NULL){
if(class(N1) != "biomass_distribution"){
	stop("N1 has to be of class 'biomass_distribution'!")
}
if(class(N2) != "biomass_distribution"){
	stop("N2 has to be of class 'biomass_distribution'!")
}
times = N1$times
if(any(N2$times != times)){
	stop("The time scales of the two distributions do not agree!")
}

N1 = N1$bmdist
N2 = N2$bmdist

n = dim(N1)[2]
num_poly = dim(N1)[1]
if(any(dim(N1) != dim(N2))) warning("biomass_compare recieved series of different dimensions!")

if(graphics == TRUE){
	plot(times,N1[1,],"l",lwd=1,lty=1,ylim=c(min(N1,N2),max(N1,N2)),main="Comparison of biomass distribution",xlab="Sample time points (days)",ylab="Amount of biomass (% of overall biomass)")
	lines(times,N2[1,],lty=1,"l",col=8,lwd=1)
	for (i in 2:num_poly){
		lines(times,N1[i,],lty=i,lwd=1)
		lines(times,N2[i,],lty=i,col=8,"l",lwd=1)
	}
	if(! is.null(legend)) legend(legend,paste("polygon",1:num_poly),fill=1:num_poly)
}


# Return relative difference
return(c(sum((N1 - N2)^2)/n,sum(abs(N1 - N2)/ifelse((N1+N2) == 0, 1, N1+N2)) / n))
}