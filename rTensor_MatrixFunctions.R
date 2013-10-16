#####(rTensor) Tensor Algebra and Statistical Models####
#####Functions that operate on Matrices and Arrays
###Hamadard (element-wise) product of a list of matrices
hamadard_list <- function(L){
	isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
	stopifnot(all(unlist(lapply(L,isvecORmat))))
	retmat <- L[[1]]
	for (i in 2:length(L)){
		retmat <- retmat*L[[i]]
	}
	retmat
}
###Kronecker product of a list of matrices
kronecker_list <- function(L){
	isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
	stopifnot(all(unlist(lapply(L,isvecORmat))))
	retmat <- L[[1]]
	for(i in 2:length(L)){
		retmat <- kronecker(retmat,L[[i]])
	}
	retmat
}
###Khatri Rao product of matrices
khatri_rao <- function(x,y){
	if (!(is.matrix(x)&&is.matrix(y))) stop("Arguments must be matrices.")
	if (dim(x)[2]!=dim(y)[2]) stop("Arguments must have same number of columns.")
	retmat <- matrix(0,nrow=dim(x)[1]*dim(y)[1],ncol=dim(x)[2])
	for (j in 1:ncol(retmat)) retmat[,j] <- kronecker(x[,j],y[,j])
	retmat
}
###Khatri Rao product of a list of matrices
khatri_rao_list <- function(L,reverse=FALSE){
	stopifnot(all(unlist(lapply(L,is.matrix))))
	ncols <- unlist(lapply(L,ncol))
	stopifnot(length(unique(ncols))==1)
	ncols <- ncols[1]
	nrows <- unlist(lapply(L,nrow))
	retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
	
	if (reverse) L <- rev(L)
	for(j in 1:ncols){
			Lj <- lapply(L,function(x) x[,j])
			retmat[,j] <- kronecker_list(Lj)
	}
	retmat
}
###Norm of vectors
norm_vec <- function(vec){
	norm(as.matrix(vec))
}
###Circulant Matrix from a vector
circ_mat <- function(vec){
	stopifnot(is.vector(vec))
	n <- length(vec)
	suppressWarnings(matrix(vec[t(matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n])],n,n))
}