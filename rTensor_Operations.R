###Tensor times matrix (m-mode product)
ttm<-function(tnsr,mat,m=NULL){
	stopifnot(is.matrix(mat))
	if(is.null(m)) stop("m must be specified")
	mat_dims <- dim(mat)
	modes_in <- getModes(tnsr)
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)
	retarr_m <- mat%*%tnsr_m
	rs_fold(retarr_m,m=m,modes=modes_out)
}
###Tensor times vector (contracted m-mode product)
ttv<-function(tnsr,vec,m=NULL){
	if(is.null(m)) stop("m must be specified")
	vec_dim <- length(vec)
	modes_in <- getModes(tnsr)	
	stopifnot(modes_in[m]==vec_dim)
	modes_out <- modes_in
	modes_out[m] <- 1
	tnsr_m <- rs_unfold(tnsr,m=m)
	retarr_m <- vec%*%tnsr_m
	rs_fold(retarr_m,m=m,modes=modes_out)
}
###Tensor times a list of matrices
ttl<-function(tnsr,list_mat,ms=NULL){
	if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
	if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
	num_mats <- length(list_mat)
	if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
	mat_nrows <- vector("list", num_mats)
	mat_ncols <- vector("list", num_mats)
	for(i in 1:num_mats){
	mat <- list_mat[[i]]
	m <- ms[i]
	mat_dims <- dim(mat)
	modes_in <- getModes(tnsr)
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)
	retarr_m <- mat%*%tnsr_m
	tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
	}	
	tnsr	
}
###Tensor Transpose
setMethod("t",signature="Tensor",
definition=function(x){
	if(getNumModes(x)!=3) stop("Tensor Transpose currently only implemented for 3d Tensors")
	modes <- getModes(x)
	new_arr <- array(apply(getData(x)[,,c(1L,modes[3]:2L)],MARGIN=3,FUN=t),dim=modes)
	as.tensor(new_arr)
})
#ifft function definition
ifft <- function(x){as.numeric(fft(x,inverse=TRUE)/length(x))}
###Tensor Multiplication (only defined for 3-d so far)
setMethod("%*%",signature=c("Tensor","Tensor"),
definition=function(x,y){
	tensor_product3d(x,y)
})
#
tensor_product3d <- function(x,y){
	if((getNumModes(x)!=3)||(getNumModes(y)!=3)) stop("Tensor Multiplication currently only implemented for 3d Tensors")
	modes_x <- getModes(x)
	modes_y <- getModes(y)
	if(modes_x[2]!=modes_y[1]) stop("Mode 2 of x and Mode 1 of y must match")
	n3 <- modes_x[3]
	if(n3!=modes_y[3]) stop("Modes 3 of x and y must match")
	#fft's for x and y
	fft_x <- aperm(apply(getData(x),MARGIN=1:2,fft),c(2,3,1))
	fft_y <- aperm(apply(getData(y),MARGIN=1:2,fft),c(2,3,1))
	#multiply the faces (this is terribad! TO-DO: think of better way!)
	fft_ret <- array(0,dim=c(modes_x[1],modes_y[2],n3))
	for(i in 1:n3){
		fft_ret[,,i]<-fft_x[,,i]%*%fft_y[,,i]
	}
	#ifft and return as Tensor
	as.tensor(aperm(apply(fft_ret,MARGIN=1:2,ifft),c(2,3,1)))
}


