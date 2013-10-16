#####(rTensor) Tensor Algebra and Statistical Models####
#####Method Definitions
#modes getter
setMethod(f="getModes",
signature="Tensor",
definition=function(x){
	x@modes
})
#num_modes getter
setMethod(f="getNumModes",
signature="Tensor",
definition=function(x){
	x@num_modes
})
#modenames getter
setMethod(f="getModenames",
signature="Tensor",
definition=function(x){
	x@modenames
})
#data getter
setMethod(f="getData",
signature="Tensor",
definition=function(x){
	if(x@num_modes==1) return(as.vector(x@data))
	return(x@data)
})
###Show and Print
setMethod(f="show",
signature="Tensor",
definition=function(x){
	cat("Numeric Tensor of", x@num_modes, "Modes\n", sep=" ")
	cat("Modes: ", x@modes, "\n", sep=" ")
	modenames <-x@modenames
	if (is.null(modenames)){
		cat("Modenames: <empty>\n")
	}else{
		cat("Modenames: ", modenames, "\n", sep=" ")
	}
	cat("Data: \n")
	print(head(x@data))
})
#
setMethod(f="print",
signature="Tensor",
definition=function(x,...){
	show(x)
})
###Head and Tail
setMethod(f="head",
signature="Tensor",
definition=function(x,...){
	head(x@data,...)
})
#
setMethod(f="tail",
signature="Tensor",
definition=function(x,...){
	tail(x@data,...)
})
###Ops
setMethod("Ops", signature(e1="Tensor", e2="Tensor"),
definition=function(e1,e2){
	e1@data<-callGeneric(e1@data, e2@data)
	validObject(e1)
	e1
})
#
setMethod("Ops", signature(e1="Tensor", e2="array"),
definition=function(e1,e2){
	e1@data<-callGeneric(e1@data,e2)
	validObject(e1)
	e1
})
#
setMethod("Ops", signature(e1="array", e2="Tensor"),
definition=function(e1,e2){
	e2@data<-callGeneric(e1,e2@data)
	validObject(e2)
	e2
})
#
setMethod("Ops", signature(e1="Tensor", e2="numeric"),
definition=function(e1,e2){
	e1@data<-callGeneric(e1@data,e2)
	validObject(e1)
	e1
})
#
setMethod("Ops", signature(e1="numeric", e2="Tensor"),
definition=function(e1,e2){
	e2@data<-callGeneric(e1,e2@data)
	validObject(e2)
	e2
})
#####Subsetting Methods ('[' defined for the array data)
setMethod("[", signature="Tensor",
definition=function(x,i,j,...){
	as.tensor(`[`(x@data,i,j,...))
})
###Sum/mean aross a given mode
setMethod("modeSum",signature="Tensor",
definition=function(x,m=NULL){
	if(is.null(m)) stop("must specify mode m")
	num_modes <- x@num_modes
	if(m<1||m>num_modes) stop("m out of bounds")
	perm <- c(m,(1L:num_modes)[-m])
	arr <- colSums(aperm(x@data,perm),dims=1L)
	as.tensor(arr)
})
#
setMethod("modeMean",signature="Tensor",
definition=function(x,m=NULL){
	if(is.null(m)) stop("must specify mode m")
	num_modes <- x@num_modes
	if(m<1||m>num_modes) stop("m out of bounds")
	perm <- c(m,(1L:num_modes)[-m])
	arr <- colSums(aperm(x@data,perm),dims=1L)
	modes <- x@modes
	as.tensor(arr/modes[m])
})
###Sweep
setMethod("sweep", signature="Tensor",
definition=function(x,m=NULL,stats=NULL,func=NULL,...){
	if(is.null(m)) stop("must specify mode m")
	as.tensor(sweep(x@data,MARGIN=m,STATS=stats,FUN=func,...))
})
###Norm and Inner Product
setMethod("fnorm",signature="Tensor",
definition=function(x){
	arr<-x@data
	sqrt(sum(arr*arr))
})
#
setMethod("innerProd",signature=c(x1="Tensor", x2="Tensor"),
definition=function(x1,x2){
	stopifnot(x1@modes==x2@modes)
	arr1 <- x1@data
	arr2 <- x2@data
	sum(as.numeric(arr1*arr2))
})























############OLD;REQUIRES ABIND##################
# setMethod("getSubtensor", signature="Tensor",
# definition=function(x,indices=NULL,dims=seq(len=max(getNumModes(x),1)),drop=NULL){
# stopifnot(require(abind))
# if(!is.vector(indices)) stop("indices must be a vector")
# idx <- as.list(indices)
# null_ind <- which(idx < 0)
# if(length(null_ind)!=0) idx[null_ind] <- list(NULL)
# subTensor <- tensor(abind::asub(getData(x),idx,dims=dims,drop=drop),modenames=getModenames(x))
# subTensor
# })
# #
# setMethod("getFiber", signature="Tensor",
# definition=function(x,indices=NULL,asTensor=FALSE){
# stopifnot(require(abind))
# if(!is.vector(indices)) stop("indices must be a vector")
# num_modes <- getNumModes(x)
# if(num_modes==1){
	# if(asTensor) return(x)
	# return(getData(x))
	# }
# if(length(indices)!=num_modes) stop ("indices must have length N")
# idx <- as.list(indices)
# null_ind <- which(idx < 0)
# if(length(null_ind)!=1) stop("there must be exactly 1 negative index")
# idx[null_ind] <- list(NULL)
# fiber <- abind::asub(getData(x),idx=idx,drop=TRUE)
# if(asTensor) fiber <- as.tensor(fiber)
# fiber
# })
# #
# setMethod("getSlice", signature="Tensor",
# definition=function(x,indices=NULL,asTensor=FALSE){
# stopifnot(require(abind))
# if(!is.vector(indices)) stop("indices must be a vector")
# num_modes <- getNumModes(x)
# if(num_modes==1) stop("Tensor has only 1 mode")
# if(num_modes==2){
	# if(asTensor) return(x)	
	# return(getData(x))		
# }
# if(length(indices)!=num_modes) stop ("indices must have length N")
# idx <- as.list(indices)
# null_ind <- which(idx < 0)
# if(length(null_ind)!=2) stop("there must be exactly 2 negative indice")
# idx[null_ind] <- list(NULL)
# slice <- abind::asub(getData(x),idx=idx,drop=TRUE)
# if (asTensor) slice <- as.tensor(slice)
# slice
# })