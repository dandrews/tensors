#####(rTensor) Tensor Algebra and Statistical Models####
#####Class Defintions
### Master Class (Made virtual so that dense tensor & sparse tensor implementations could easily inherit basic tensor properties)
setClass("Tensor",
representation(num_modes = "integer", modes = "integer", modenames = "list", "VIRTUAL"),
validity = function(object){
	num_modes <- object@num_modes
	modes <- object@modes
	modenames <- object@modenames
	errors <- character()
	if ((length(num_modes)==0)||(length(modes)==0)) return(TRUE) #empty object
	if (any(modes <= 0)){
	msg <- "'modes' must contain positive values; if any mode is 0, consider a smaller num_modes"
	errors <- c(errors, msg)
	}
	if (!is.list(modenames)){
		msg <- "'modenames' must be a list"
	}
	if (!is.null(modenames[[1]]) && (length(modenames)!=num_modes)){
	msg <- "warning: 'modenames' length does not match number of modes. recycling"
	cat(msg)
	}
	if(length(errors)==0) TRUE else errors
})
###Dense Tensor Classes (Non-virtual)
### Numeric Dense Tensor
setClass("ndTensor", representation(data = "array"), contains = "Tensor", 
#prototype = prototype(data = array(data=1)),
validity = function(object){
	errors <- character()
	ndtv <- "Pass" #ndtv <- ndTensor_validate(object@data)
	if(ndtv!="Pass"){
	errors <- c(errors, ndtv)
	}	
})
#####Initializations 
setMethod(f="initialize",
signature="ndTensor",
definition = function(.Object, num_modes=NULL, modes=NULL, modenames=list(NULL), data=NULL){
	if(is.null(data)){
		.Object@num_modes <- integer(0)
		.Object@modes <- integer(0)
		.Object@modenames <- list(NULL)
		validObject(.Object)
		.Object
	}
	if(is.null(num_modes)){
		if (is.vector(data)) num_modes <- 1L
		else{num_modes <- length(dim(data))}
	}
	if(is.null(modes)){
		if (is.vector(data)) modes <- length(data)
		else{modes <- dim(data) }
	}
	if(is.null(modenames[[1]])&&!is.null(dimnames(.Object))){
		modenames <- dimnames(data)
		}
	.Object@num_modes <- num_modes
	.Object@modes <- modes
	.Object@modenames <- modenames
	.Object@data <- as.array(data,dim=modes,dimnames=modenames)
	validObject(.Object)
	.Object
})


















##################################<to-do-list>######################################
#implement tensor normal!!!
#########################<to be implemented MUCH later>#############################
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
####################################################################################