#####(rTensor) Tensor Algebra and Statistical Models####
#####Class Defintions (Only 1 class in version 0.1)
setClass("Tensor",
representation(num_modes = "integer", modes = "integer", modenames = "ANY", data="array"),
validity = function(object){
	num_modes <- object@num_modes
	modes <- object@modes
	modenames <- object@modenames
	errors <- character()
	if (any(modes <= 0)){
		msg <- "'modes' must contain strictly positive values; if any mode is 1, consider a smaller num_modes"
		errors <- c(errors, msg)
	}
	if (!is.null(modenames) && (length(modenames)!=num_modes)){
		msg <- "warning: 'modenames' length does not match number of modes. recycling"
		errors <- c(errors, msg)
	}
	if(length(errors)==0) TRUE else errors
})

#####Initialization
setMethod(f="initialize",
signature="Tensor",
definition = function(.Object, num_modes=NULL, modes=NULL, modenames=NULL, data=NULL){
	if(is.null(num_modes)){
		if (is.vector(data)) num_modes <- 1L
		else{num_modes <- length(dim(data))}
	}
	if(is.null(modes)){
		if (is.vector(data)) modes <- length(data)
		else{modes <- dim(data)}
	}
	if(is.null(modenames)&&!is.null(dimnames(.Object))){
		modenames <- dimnames(data)
		}
	.Object@num_modes <- num_modes
	.Object@modes <- modes
	.Object@modenames <- modenames
	.Object@data <- array(data,dim=modes,dimnames=modenames)
	validObject(.Object)
	.Object
})

#####Creation of ndTensor from an array/matrix/vector
as.tensor <- function(x, mode=NULL, modenames=NULL,drop=TRUE){
	stopifnot(is.array(x)||is.vector(x))
	if (is.vector(x)){
		modes <- c(length(x))
		num_modes <- 1L
	}else{
		modes <- dim(x)
		num_modes <- length(modes)
		if(is.null(modenames)&&!is.null(dimnames(x))){
			modenames <- dimnames(data)
			}
		dim1s <- which(modes==1)
		if(drop && (length(dim1s)>0)){
			modes <- modes[-dim1s]
			if(!is.null(modenames[[1]])) modenames <- modenames[-dim1s]
			num_modes <- num_modes-length(dim1s)
		}
	}
new("Tensor",num_modes,modes,modenames,data=array(x,dim=modes,dimnames=modenames))
}
