tensors
=======

R package for multi-linear statistics.

A collaboration between James Li and Daniel Andrews

##################################<to-do-list>######################################
#########################<to be implemented MUCH later>#############################
#implement tensor normal
#(Make Tensor virtual so that dense tensor & sparse tensor implementations could easily inherit basic tensor properties)
#Example: Dense Tensor Classes (Non-virtual)
### Numeric Dense Tensor
#setClass("ndTensor", representation(data = "array"), contains = "Tensor" 
#validity = function(object){
#	errors <- character()
#	ndtv <- "Pass" #ndtv <- ndTensor_validate(object@data)
#	if(ndtv!="Pass"){
#	errors <- c(errors, ndtv)
#	}	
#})
### Integer Dense Tensor
#setClass("idTensor", representation(data = "integer"), contains = c("Tensor", "array"), validity = function(object) idTensor_validate(object))
### Logical Dense Tensor
#setClass("ldTensor", representation(data = "logical"), contains = c("Tensor", "array"), validity = function(object) ldTensor_validate(object))
### Sparse tensor

### Initialization Functions
#setMethod("initialize", "ndTensor", function(.Object){
#	modes <- 
#	.Object@data <- .Internal(array(data,modes,modenames)))
#	}

### Rcpp & marry.hxx 
#require(Rcpp)
#Sys.setenv("PKG_CXXFLAGS"="-I /Users/jamesyili/cpp_include/")
#sourceCpp(file="/Users/jamesyili/Dropbox/Advanced R/tensor.cpp")
#
#a <- array(1:32, dim=rep(2,5))
#b <- marrayC(a, dim(a))
#
#dims = c(1,24,5,12,12,9,8)
#size = prod(dims)
#a <- array(runif(size), dim = dims)

###Subset setters (NEED refClass??)
# setMethod("setSubtensor", signature="ndTensor",
# definition=function(x, dim=seq(len=max(getNumModes(x),1))){
# })
# setMethod("setFiber", signature="ndTensor",
# definition=function(x, dim=seq(len=max(getNumModes(x),1))){
# })
# setMethod("setSlice", signature="ndTensor",
# definition=function(x, dim=seq(len=max(getNumModes(x),1))){
# })

####################################################################################