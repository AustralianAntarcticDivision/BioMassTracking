#' BiomassTracking
#' 
#' @name BiomassTracking
#' @docType package
#' @useDynLib BiomassTracking, .registration = TRUE
#' @importFrom grDevices dev.new nclass.FD nclass.Sturges nclass.scott terrain.colors topo.colors 
#' @importFrom graphics arrows filled.contour image lines matplot par plot points
#' @importFrom stats na.omit quantile rnorm runif
#' @importFrom utils read.table
NULL

#' Polygon structure and land positions
#' 
#' Example data set of a polygon structure over the Kerguelen plateau. 0
#' denotes land and positive integer values give the number of the polygon the
#' respective grid point belongs to.
#' 
#' 
#' @name Sdata
#' @docType data
#' @format A matrix of integer values.
#' @source Thorsten Lenser, just an example.
#' @keywords datasets
NULL





#' EW component of a flow field
#' 
#' Example data set of u-components of the flow field over the Kerguelen
#' plateau. Derived from satellite measurements of sea surface height.
#' 
#' 
#' @name Udata
#' @docType data
#' @format Three columns, containing longitude, latitude and current speed in
#' m/s. Unknown values are given as -9999.
#' @source Belinda Ronai, Australian Antarctic Division, unpublished data.
#' @keywords datasets
NULL





#' NS component of a flow field
#' 
#' Example data set of v-components of the flow field over the Kerguelen
#' plateau. Derived from satellite measurements of sea surface height.
#' 
#' 
#' @name Vdata
#' @docType data
#' @format Three columns, containing longitude, latitude and current speed in
#' m/s. Unknown values are given as -9999.
#' @source Belinda Ronai, Australian Antarctic Division, unpublished data.
#' @keywords datasets
NULL



