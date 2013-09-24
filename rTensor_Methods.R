#####(rTensor) Tensor Algebra and Statistical Models####

#####Method Definitions

###Accessor (Getter)
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
	if(getNumModes(x)==1) return(as.vector(x@data))
	return(x@data)
})

###Show and Print
setMethod(f="show",
signature="Tensor",
definition=function(x){
	cat("Numeric Tensor of", getNumModes(x), "Modes\n", sep=" ")
	cat("Modes: ", getModes(x), "\n", sep=" ")
	modenames <- getModenames(x)
	if (is.null(modenames[[1]])){
		cat("Modenames: <empty>\n")
	}else{
		cat("Modenames: ", unlist(modenames), "\n", sep=" ")
	}
	cat("Data: \n")
	print(head(getData(x)))
})

setMethod(f="print",
signature="Tensor",
definition=function(x,...){
	show(x)
})

###Head and Tail
setMethod(f="head",
signature="Tensor",
definition=function(x,...){
	show(x)
	invisible(head(getData(x)))
})

setMethod(f="head",
signature="Tensor",
definition=function(x,...){
	show(x)
	invisible(tail(getData(x)))
})

###Ops
setMethod("Ops", signature(e1="Tensor", e2="Tensor"),
function(e1, e2){
	e1@data<-callGeneric(getData(e1), getData(e2))
	validObject(e1)
	e1
})

setMethod("Ops", signature(e1="Tensor", e2="array"),
function(e1,e2){
	e1@data<-callGeneric(geData(e1),e2)
	validObject(e1)
	e1
})

setMethod("Ops", signature(e1="array", e2="Tensor"),
function(e1,e2){
	e2@data<-callGeneric(e1,getData(e2))
	validObject(e2)
	e2
})

setMethod("Ops", signature(e1="Tensor", e2="numeric"),
function(e1,e2){
	e1@data<-callGeneric(getData(e1),e2)
	validObject(e1)
	e1
})

setMethod("Ops", signature(e1="numeric", e2="Tensor"),
function(e1,e2){
	e2@data<-callGeneric(e1,getData(e2))
	validObject(e2)
	e2
})

#####Subsetting Methods
###'[' and '[<-' defined for the array data
setMethod("[", signature="Tensor",
definition=function(x,...,drop=FALSE){
	'['(getData(x),...,drop=drop)
})

setReplaceMethod("[",signature="Tensor",
definition=function(x,...,value){
	as.tensor('[<-'(getData(x),...,value),modenames=getModenames(x))
})

###Subset getters
setMethod("getSubtensor", signature="Tensor",
definition=function(x,indices=NULL,dims=seq(len=max(getNumModes(x),1)),drop=NULL){
	stopifnot(require(abind))
	if(!is.vector(indices)) stop("indices must be a vector")
	idx <- as.list(indices)
	null_ind <- which(idx < 0)
	if(length(null_ind)!=0) idx[null_ind] <- list(NULL)
	subTensor <- tensor(abind::asub(getData(x),idx,dims=dims,drop=drop),modenames=getModenames(x))
	subTensor
})

setMethod("getFiber", signature="Tensor",
definition=function(x,indices=NULL,asTensor=FALSE){
	stopifnot(require(abind))
	if(!is.vector(indices)) stop("indices must be a vector")
	num_modes <- getNumModes(x)
	if(num_modes==1){
		if(asTensor) return(x)
		return(getData(x))
		}
	if(length(indices)!=num_modes) stop ("indices must have length N")
	idx <- as.list(indices)
	null_ind <- which(idx < 0)
	if(length(null_ind)!=1) stop("there must be exactly 1 negative index")
	idx[null_ind] <- list(NULL)
	fiber <- abind::asub(getData(x),idx=idx,drop=TRUE)
	if(asTensor) fiber <- as.tensor(fiber)
	fiber
})

setMethod("getSlice", signature="Tensor",
definition=function(x,indices=NULL,asTensor=FALSE){
	stopifnot(require(abind))
	if(!is.vector(indices)) stop("indices must be a vector")
	num_modes <- getNumModes(x)
	if(num_modes==1) stop("Tensor has only 1 mode")
	if(num_modes==2){
		if(asTensor) return(x)	
		return(getData(x))		
	}
	if(length(indices)!=num_modes) stop ("indices must have length N")
	idx <- as.list(indices)
	null_ind <- which(idx < 0)
	if(length(null_ind)!=2) stop("there must be exactly 2 negative indice")
	idx[null_ind] <- list(NULL)
	slice <- abind::asub(getData(x),idx=idx,drop=TRUE)
	if (asTensor) slice <- as.tensor(slice)
	slice
})

#####Tensor Unfoldings
###Matricization (unfolding) in the m mode - aka Row Space Unfolding
###NEW DEFN
setMethod("rs_unfold", signature="Tensor",
definition=function(x,m=NULL,asTensor=FALSE){
	rs <- m
	cs <- (1:num_modes)[-m]
	g_unfold(x,rs=rs,cs=cs,asTensor=asTensor)
})
###OLD DEFN
# setMethod("m_unfold", signature="Tensor",
# definition=function(x,m=NULL,asTensor=FALSE){
	# if(is.null(m)) stop("mode m must be specified")
	# num_modes <- getNumModes(x)
	# if(m < 1 || m > num_modes) stop("mode m incorrectly specified")
	# modes <- getModes(x)
	# mat <- getData(x)
	# new_modes <- c(modes[m],prod(modes[-m]))
	# if(m == 1) {
		# dim(mat) <- new_modes
		# if(asTensor) mat <- as.tensor(mat)
		# return(mat)
	# }
	# if(m == num_modes){
		# perm <- c(m,1:(m-1))
	# }else {
		# perm <- c(m,1:(m-1),(m+1):num_modes)
	# }
	# mat <- aperm(mat,perm)
	# dim(mat) <- new_modes
	# if(asTensor) mat <- as.tensor(mat)
	# mat
# })
###See rTensor_Functions for Un-matricization (not a method since it operates on matrices, not tensors)

###Column Space Unfolding
setMethod("cs_unfold", signature="Tensor",
definition=function(x,m=NULL,asTensor=FALSE){
	num_modes <- getNumModes(x)
	rs <- (1:num_modes)[-m]
	cs <- m
	unfold(x,rs=rs,cs=cs,asTensor=asTensor)
})

###General Unfolding
setMethod("unfold", signature="Tensor",
definition=function(x,rs=NULL,cs=NULL,asTensor=FALSE){
	if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
	num_modes <- getNumModes(x)
	if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
	if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
	perm <- c(rs,cs)
	if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
	modes <- getModes(x)
	mat <- getData(x)
	new_modes <- c(prod(modes[rs]),prod(modes[cs]))
	mat <- aperm(mat,perm)
	dim(mat) <- new_modes
	if(asTensor) mat <- as.tensor(mat)
	mat
})
###See rTensor_Functions for Un-matricization (not a method since it operates on matrices, not tensors)

###Sweep
setMethod("sweep", signature="Tensor",
definition=function(x,m=NULL,operand=NULL,operator=NULL,asTensor=TRUE,...){
	if(is.null(m)) stop("must specify mode m")
	arr <- getData(x)
	modenames <- getModenames(x)
	arr <- sweep(arr,MARGIN=m,STATS=stats,FUN=func,...)
	if(asTensor) return(as.tensor(arr,modenames=modenames))
	arr
})


###Norm and Inner Product
setMethod("fnorm",signature="Tensor",
definition=function(x){
	arr<-getData(x)
	sqrt(sum(arr*arr))
})

setMethod("inner_prod",signature=c(e1="Tensor", e2="Tensor"),
definition=function(e1,e2){
	stopifnot(getModes(e1)==getModes(e2))
	arr1 <- getData(e1)
	arr2 <- getData(e2)
	sum(arr1*arr2)
})

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





