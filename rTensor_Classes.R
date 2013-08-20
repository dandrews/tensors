#####(rTensor) Tensor Algebra and Statistical Models####

#####Class Defintions
### Master Class (Made virtual so that dense tensor & sparse tensor implementations could easily inherit basic tensor properties)
setClass("Tensor",
representation(num_modes = "integer", modes = "integer", modenames = "list", "VIRTUAL"),
#prototype = prototype(num_modes = 1L, modes = c(1L), modenames= list(NULL)),
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
	if (is.null(modenames)){modenames <- list(NULL)
		}else if (!is.list(modenames)){
		msg <- "if non-empty, 'modenames' must be a list"
	}
	if (is.list(modenames) && !is.null(modenames[[1]]) && (length(modenames)!=num_modes)){
	msg <- "warning: 'modenames' length does not match number of modes; recycling."
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
definition = function(.Object, num_modes=NULL, modes=NULL, modenames=NULL, data=NULL){
cat("~~~Initializing a new ndTensor~~~~ \n")
	if(is.null(data)){
		.Object@num_modes <- integer(0)
		.Object@modes <- integer(0)
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
	if(!is.null(modenames)){
		modenames <- dimnames(data)
	}else{ modenames <- list(NULL)}
	
	.Object@num_modes <- num_modes
	.Object@modes <- modes
	.Object@modenames <- modenames
	.Object@data <- array(data)
	validObject(.Object)
	.Object
})

###Helper Function to Create ndTensor
as.tensor <- function(x, num_modes, modes, modenames=NULL, type = "numeric"){
if(!(type %in% c("numeric", "integer", "logical"))) stop("type must be 'numeric', 'integer', or 'logical'")
if (is.vector(x)){
	modes <- c(length(x))
	num_modes <- 1L
	modenames <- names(x)
}else if (is.matrix(x)){
	modes <- dim(x)
	num_modes <- 2L
	modenames <- dimnames(x)
}else if (is.array(x)){
	modes <- dim(x)
	num_modes <- length(modes)
	modenames <- dimnames(x)
}else{
stop("can only create a tensor from vectors, matrices, and arrays")	
}	

tnsr<-switch(type,
numeric=new("ndTensor",num_modes,modes,modenames,data=x),
integer=new("ndTensor",num_modes,modes,modenames,data=x),
logical=new("ndTensor",num_modes,modes,modenames,data=x)
)
return(tnsr)
}

#arr <- array(runif(1e3),dim=c(10,10,10))
#mat <- matrix(runif(100),nrow=10)
#vec <- 1:10
#arrT <- as.tensor(arr)
#matT <- as.tensor(mat)
#vecT <- as.tensor(vec)


###Plot Functions



##################################<to-do-list>######################################
#### in R:
#implement add/subtract
#implement tensor times vector
#implement tensor times matrix
#implement matricization
#implement un-matricization
#implement alternating least squares



#### in C++:
#implement ndTensor_validate















##########################<to be implemented MUCH later>#############################
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

