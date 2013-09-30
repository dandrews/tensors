#####(rTensor) Tensor Algebra and Statistical Models####
#####Tensor Unfoldings
###General Unfolding
setMethod("unfold", signature="Tensor",
definition=function(x,rs=NULL,cs=NULL,asTensor=FALSE){
	#checks
	if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
	num_modes <- getNumModes(x)
	if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
	if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
	perm <- c(rs,cs)
	if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
	modes <- getModes(x)
	mat <- getData(x)
	new_modes <- c(prod(modes[rs]),prod(modes[cs]))
	#rearranges into a matrix
	mat <- aperm(mat,perm)
	dim(mat) <- new_modes
	if(asTensor) mat <- as.tensor(mat)
	mat
})
###Row Space Unfolding in the m mode - aka Matricization (Kolda et. al)
setMethod("rs_unfold", signature="Tensor",
definition=function(x,m=NULL,asTensor=FALSE){
	if(is.null(m)) stop("mode m must be specified")
	num_modes <- getNumModes(x)
	rs <- m
	cs <- (1:num_modes)[-m]
	unfold(x,rs=rs,cs=cs,asTensor=asTensor)
})
###Column Space Unfolding in the m mode (Martin et. al)
setMethod("cs_unfold", signature="Tensor",
definition=function(x,m=NULL,asTensor=FALSE){
	if(is.null(m)) stop("mode m must be specified")
	num_modes <- getNumModes(x)
	rs <- (1:num_modes)[-m]
	cs <- m
	unfold(x,rs=rs,cs=cs,asTensor=asTensor)
})
#####Matrix Foldings
###General folding (inverse function to unfold)
fold <- function(mat, rs = NULL, cs = NULL, modes=NULL){
	#checks
	if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	if(!is(mat,"Tensor")){
		if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
		}else{
			stopifnot(getNumModes(mat)==2)
			mat <- getData(mat)			
			}
	num_modes <- length(modes)
	stopifnot(num_modes==length(rs)+length(rs))
	mat_modes <- dim(mat)
	if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
	#rearranges into array
	iperm <- match(1:num_modes,c(rs,cs))
	arr<-array(mat,dim=c(modes[rs],modes[cs]))
	as.tensor(aperm(arr,iperm))
}
###Row Space Folding (inverse funtion to rs_unfold) in the m mode
rs_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	rs <- m
	cs <- (1:num_modes)[-m]
	fold(mat,rs=rs,cs=cs,modes=modes)
}
###Col Space Folding (inverse function to cs_unfold) in the m mode
cs_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	cs <- m
	rs <- (1:num_modes)[-m]
	fold(mat,rs=rs,cs=cs,modes=modes)	
}


















###OLD DEFN of rs_unfold: DO NOT USE###
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