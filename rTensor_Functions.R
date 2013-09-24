#####(rTensor) Tensor Algebra and Statistical Models####

#####Tensor Functions
###Creation of ndTensor
tensor <- function(x, modes, modenames=list(NULL), type = "numeric",drop=TRUE){
	as.tensor(x=x, modes=modes, modenames=modenames, type = type,drop=drop)
}

as.tensor <- function(x, modes, modenames=list(NULL), type = "numeric",drop=TRUE){
if(!(type %in% c("numeric", "integer", "logical"))) stop("type must be 'numeric', 'integer', or 'logical'")
if (is.vector(x)){
	modes <- c(length(x))
	num_modes <- 1L
}else if (is.array(x)){
	modes <- dim(x)
	num_modes <- length(modes)
	if(is.null(modenames[[1]])&&!is.null(dimnames(x))){
		modenames <- dimnames(data)
		}
	dim1s <- which(modes==1)
	if(drop && (length(dim1s)>0)){
		modes <- modes[-dim1s]
		if(!is.null(modenames[[1]])) modenames <- modenames[-dim1s]
		num_modes <- num_modes-length(dim1s)
	}
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

###Un-matricization for rs_unfold
rs_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	if(!is(mat,"Tensor")){
		if(!(class(mat) %in% c("array", "matrix")))  stop("mat must be of class 'array', 'matrix'")
		}else{
			mat <- getData(mat)			
			}
	if((class(mat)=="array") && (dim(mat) != 2)) stop("mat must have 2 dimensions")
	num_modes <- length(modes)
	# if(num_modes < 3) stop("Tensor must have 3 or more modes")
	if(m < 1 || m > num_modes) stop("mode m incorrectly specified")
	mat_modes <- dim(mat)
	if((mat_modes[1]!=modes[m]) || (mat_modes[2]!=prod(modes[-m]))) stop("matrix nrow/ncol does not match Tensor modes")
	if(m==1){
		return(as.tensor(array(mat,dim=modes)))
	}
	new_modes <- c(modes[m],modes[-m])
	arr <- array(mat,dim=new_modes)
	if (m == num_modes){
		iperm <- c(2:num_modes,1)
	}else{
		iperm <- c(2:m,1,(m+1):num_modes)
		}
	as.tensor(aperm(arr,iperm))
}

cs_fold <- function(){
	
}

###General folding
fold <- function(mat, rs = NULL, cs = NULL){
	
}

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




###Tensor times matrix (m-mode product)
ttm <- function(tnsr, mat, m=NULL){
	stopifnot(is(tnsr,"Tensor") && is.matrix(mat))
	if(is.null(m)) stop("m must be specified")
	mat_dims <- dim(mat)
	modes_in <- getModes(tnsr)
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- m_unfold(tnsr,m=m)
	retarr_m <- mat%*%tnsr_m
	m_fold(retarr_m,m=m,modes=modes_out)
}

###Tensor times vector (contracted m-mode product)
ttv <- function(tnsr,vec,m=NULL){
	stopifnot(is(tnsr,"Tensor") && is.vector(vec))
	if(is.null(m)) stop("m must be specified")
	vec_dim <- length(vec)
	modes_in <- getModes(tnsr)	
	stopifnot(modes_in[m]==vec_dim)
	modes_out <- modes_in
	modes_out[m] <- 1
	tnsr_m <- m_unfold(tnsr,m=m)
	retarr_m <- vec%*%tnsr_m
	m_fold(retarr_m,m=m,modes=modes_out)
}

###Tensor times a list of matrices
ttl <- function(tnsr, list_mat, ms = NULL){
	stopifnot(is(tnsr,"Tensor"))
	stopifnot(is.list(list_mat))
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
	tnsr_m <- m_unfold(tnsr,m=m)
	retarr_m <- mat%*%tnsr_m
	tnsr <- m_fold(retarr_m,m=m,modes=modes_out)
	}
	
	tnsr
}

###Norm of vectors
norm_vec <- function(vec){
	norm(as.matrix(vec))
}

#################
# "tensor" <-
# function(A, B, alongA = integer(0), alongB = integer(0))
# {
  # A <- as.array(A)
  # dimA <- dim(A)
  # dnA <- dimnames(A)
  # if (nnA <- is.null(dnA))
    # dnA <- rep(list(NULL), length(dimA))

  # B <- as.array(B)
  # dimB <- dim(B)
  # dnB <- dimnames(B)
  # if (nnB <- is.null(dnB))
    # dnB <- rep(list(NULL), length(dimB))

  # if (length(alongA) != length(alongB))
    # stop("\"along\" vectors must be same length")

  # # special case of both length zero

  # if (length(alongA) == 0) {
    # R <- as.vector(A) %*% t(as.vector(B))
    # dim(R) <- c(dimA, dimB)
    # if (!(nnA && nnB))
    # dimnames(R) <- c(dnA, dnB)
    # return(R)
  # }

  # mtch <- dimA[alongA] == dimB[alongB]
  # if (any(is.na(mtch)) || !all(mtch))
    # stop("Mismatch in \"along\" dimensions")

  # seqA <- seq(along=dimA)
  # allA <- length(seqA) == length(alongA)
  # permA <- c(seqA[-alongA], alongA)
  # if (!all(seqA == permA))
    # A <- aperm(A, permA)
  # dim(A) <- c(
    # if (allA) 1 else prod(dimA[-alongA]),
    # prod(dimA[alongA])
  # )

  # seqB <- seq(along=dimB)
  # allB <- length(seqB) == length(alongB)
  # permB <- c(alongB, seqB[-alongB])
  # if (!all(seqB == permB))
    # B <- aperm(B, permB)
  # dim(B) <- c(
    # prod(dimB[alongB]),
    # if (allB) 1 else prod(dimB[-alongB])
  # )

  # R <- A %*% B

  # if (allA && allB)
    # R <- drop(R)
  # else {
    # dim(R) <- c(
      # if (allA) integer(0) else dimA[-alongA],
      # if (allB) integer(0) else dimB[-alongB]
    # )
    # if (!(nnA && nnB))
      # dimnames(R) <- c(dnA[-alongA], dnB[-alongB])
  # }
  # R
# }

# "%*t%" <- function(x, y) tensor(x, y, 2, 2)

# "%t*%" <- function(x, y) tensor(x, y, 1, 1)

# "%t*t%" <- function(x, y) tensor(x, y, 1, 2)

