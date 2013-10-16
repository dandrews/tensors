#####(rTensor) Tensor Algebra and Statistical Models####
#####Various Tensor Decompositions

###(Truncated-)Higher-order SVD
#Input: Tensor tnsr, ranks = k_1 ,..., k_M (optional)
#Output: List containing a core tensor Z (all-orthogonal), and a list of matrices of sizes I_m x k_M if ranks provided, o/w I_m x I_m
hosvd <- function(tnsr,ranks=NULL){
	#stopifnot(is(tnsr,"Tensor"))
	num_modes <- tnsr@num_modes
	#no truncation if ranks not provided
	if(is.null(ranks)){
		cat("!ranks not provided so left singular matrices will not be truncated.\n")
		ranks <- tnsr@modes
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=num_modes,style=3)
	#loops through and performs SVD on mode-m matricization of tnsr
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		U_list[[m]] <- (svd(rs_unfold(tnsr,m=m)@data)$u)[,1:ranks[m]]
		setTxtProgressBar(pb,m)
	}
	close(pb)
	#computes the core tensor
	Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
	resid <- fnorm(ttl(Z,U_list,ms=1:num_modes)-tnsr)/fnorm(tnsr)
	#put together the return list, and returns
	list(Z=Z,U=U_list,resid=resid)	
}

###CP Decomp ALS
#Input: Tensor, number of components to decompose the tensor (k)
#Output: List containing weights lambda (vector of length k), and a list of matrices U's, were each U_m is of size I_m x k
cp_als <- function(tnsr, num_components=NULL,max_iter=500, tol=1e-6){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 1:num_modes){
		unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
	}
	est <- tnsr
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resids <- rep(0, max_iter)
	CHECK_CONV <- function(est){
		curr_resid <- fnorm(tnsr - est)
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}	
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)	
		for(m in 1:num_modes){
			V <- hamadard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
			V_inv <- solve(V)			
			tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],rev=TRUE)%*%V_inv
			lambdas <- apply(tmp,2,norm_vec)
			U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
			Z <- superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
			est <- ttl(Z,U_list,ms=1:num_modes)
		}
		#checks convergence
		if(CHECK_CONV(est)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			 }
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	invisible(list(lambdas=lambdas, U=U_list, conv=converged, norm_percent=norm_percent, resids=fnorm_resids))
}

###Tucker Decomp ALS
#Input: Tensor of size I_1,...,I_M, ranks k_1,...,k_M, where M is the number of modes to Tensor
#Output: List containing a core Tensor Z and and a list of matrices U's, were each U_m is of size I_m x k_m
tucker_als <- function(tnsr,ranks=NULL,max_iter=50,tol=1e-6){
	stopifnot(is(tnsr,"Tensor"))
	if(is.null(ranks)) stop("ranks must be specified")
	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		U_list[[m]] <- (svd(rs_unfold(tnsr,m=m)@data)$u)[,1:ranks[m]]
	}
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resids <- rep(0, max_iter)
	CHECK_CONV <- function(Z,U_list){
		est <- ttl(Z,U_list,ms=1:num_modes)
		curr_resid <- fnorm(tnsr - est)
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)	
	modes_seq <- 1:num_modes
		for(m in modes_seq){
			#core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
			#truncated SVD of X
			U_list[[m]] <- (svd(rs_unfold(X,m=m)@data)$u)[,1:ranks[m]]
		}
		#compute core tensor Z
		Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

		#checks convergence
		if(CHECK_CONV(Z, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)	
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	invisible(list(Z=Z, U=U_list, conv=converged, norm_percent = norm_percent, resids=fnorm_resids))
}

###MPCA
#Input: Tensor of size I_1,...,I_M, where m=1 is the measurement mode, ranks k_2,...,k_M, where M is the number of modes to Tensor
#Output: List containing an extended core Tensor Z_ext and and a list of matrices U's, were each U_m is of size I_m x k_m for m = 2, ..., M
mpca_als <- function(tnsr, ranks = NULL, max_iter = 500, tol=1e-6){
	if(is.null(ranks)) stop("ranks must be specified")
	stopifnot(is(tnsr,"Tensor"))
	#initialization via hosvd of M-1 modes
	num_modes <- tnsr@num_modes
	stopifnot(length(ranks)==(num_modes-1))
	ranks <- c(1,ranks)
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 2:num_modes){
		unfolded_mat <- rs_unfold(tnsr,m=m)@data
		mode_m_cov <- unfolded_mat%*%t(unfolded_mat)
		U_list[[m]] <- (svd(mode_m_cov)$u)[,1:ranks[m]]
	}
	Z_ext <- ttl(tnsr,lapply(U_list[-1],t),ms=2:num_modes)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resids <- rep(0, max_iter)
	CHECK_CONV <- function(Z_ext,U_list){
		est <- ttl(Z_ext,U_list[-1],ms=2:num_modes)
		curr_resid <- fnorm(tnsr - est)
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
	modes_seq <- 2:num_modes
		for(m in modes_seq){
			#extended core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-c(1,m)],t),ms=modes_seq[-(m-1)])
			#truncated SVD of X
			U_list[[m]] <- (svd(rs_unfold(X,m=m)@data)$u)[,1:ranks[m]]
		}
		#compute core tensor Z_ext
		Z_ext <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)
		#checks convergence
		if(CHECK_CONV(Z_ext, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	invisible(list(Z_ext=Z_ext, U=U_list, conv=converged, norm_percent = norm_percent, resids=fnorm_resids))
}

###T_SVD3d
#Input: A 3d tensor of n1 x n2 x n
#Output: List containing the U array (n1xn1xn3), V array (n2xn2xn3), S-mat (n3 x m), where m = min(n1,n2). The rows of S-mat
t_svd3d<-function(tnsr){
	if(tnsr@num_modes!=3) stop("T-SVD only implemented for 3d so far")
	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_mat <- matrix(0,nrow=m,ncol=n3)
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_mat[,j] <- decomp$d #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2,ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2,ifft),c(2,3,1)))
	S <- as.tensor(apply(S_mat,MARGIN=1,ifft))
	invisible(list(U=U,V=V,S=S))
}

###reconstruct
t_svd_reconstruct <- function(L){
	Umodes <- L$U@modes
	n1 <- Umodes[1]
	n2 <- L$V@modes[1]
	n3 <- Umodes[3]
	S_fdiagonal <- array(0,c(n1,n2,n3))
	S <- L$S@data
	for (i in 1:n3){
		S_fdiagonal[,,i] <- diag(S[i,],nrow=n1,ncol=n2)
	}
	S_fdiagonal <- as.tensor(S_fdiagonal)
#	diagcol <- function(col,n1,n2){diag(col,nrow=n1,ncol=n2)}
#	S_fdiag <- aperm(apply(S@data,MARGIN=2,diagcol,n1,n2),c(1,2,3))
	L$U%*%S_fdiagonal%*%t(L$V)
}









###T_SVD_decomp3d
#Input: A 3d tensor of n1 x n2 x n3, cutoffs k1, k2
#Output:
# t_svd_approx3d<-function(tnsr,ranks=NULL){
	# if(getNumModes(tnsr)!=3) stop("T-SVD approx only implemented for 3d so far")
	# if(is.null(ranks)||length(ranks)!=2) stop("ranks needs to be a vector of length 2")
	# modes <- tnsr@modes
	# mat <- modeSum(tnsr,m=3)@data
	# full <- svd(mat,nu=modes[1],nv=modes[2])
	# U_trunc_t <- t(full$u[,1:ranks[1]])
	# V_trunc <- full$v[,1:ranks[2]]
	# arr <- tnsr@data
	# ret_arr <- array(0,dim=c(ranks,modes[3]))
	# for (i in 1L:modes[3]){
	# ret_arr[,,i] <- U_trunc_t%*%arr[,,i]%*%V_trunc
	# }
	# as.tensor(ret_arr)
# }

###MICA - Requires W, the mixing matrix, to be specified (TO-DO)
#Input: Tensor of size I_1,...,I_M, where m=1 is the measurement mode, ranks k_2,...,k_M, where M is the number of modes to Tensor
#Output: List containing an extended core Tensor Z_ext and and a list of orthogonal matrices U's, were each U_m is of size I_m x k_m for m = 2, ..., M
# mica_als <- function(tnsr,ranks = NULL, max_iter = 500, tol=1e-6){
	# require(MASS)
	
	# if(is.null(ranks)) stop("ranks must be specified")
	# stopifnot(is(tnsr,"Tensor"))

	# #initialization via hosvd of M-1 modes
	# num_modes <- getNumModes(tnsr)
	# stopifnot(length(ranks)==(num_modes-1))
	# ranks <- c(1,ranks)
	# modes <- getModes(tnsr)
	# U_list <- vector("list",num_modes)
	# unfolded_mat <- vector("list",num_modes)
	# for(m in 2:num_modes){
		# unfolded_mat <- m_unfold(tnsr,m=m)
		# mode_m_cov <- unfolded_mat%*%t(unfolded_mat)
		# U_list[[m]] <- (svd(mode_m_cov)$u)[,1:ranks[m]]
	# }
	
	# Z_ext <- ttl(tnsr,lapply(U_list[-1],solve),ms=2:num_modes)
	# curr_iter <- 1
	# converged <- FALSE

	# #set up convergence check
	# fnorm_resids <- rep(0, max_iter)
	# CHECK_CONV <- function(Z_ext,U_list){
		# est <- ttl(Z_ext,U_list[-1],ms=2:num_modes)
		# curr_resid <- fnorm(tnsr - est)
		# cat("residual: ",curr_resid,"\n")
		# fnorm_resids[curr_iter] <<- curr_resid
		# if (curr_iter==1) return(FALSE)
		# if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		# FALSE
	# }
	
	# #main loop (until convergence or max_iter)
	# while((curr_iter < max_iter) && (!converged)){
	# cat("iteration: ",curr_iter,"\t")
	# modes_seq <- 2:num_modes
		# for(m in modes_seq){
			# #extended core Z minus mode m
			# X <- ttl(tnsr,lapply(U_list[-c(1,m)],solve),ms=modes_seq[-(m-1)])
			# #truncated SVD of X
			# U_list[[m]] <- (svd(m_unfold(X,m=m))$u)[,1:ranks[m]]
		# }
		# #compute core tensor Z_ext
		# Z_ext <- ttm(X,mat=solve(U_list[[num_modes]]),m=num_modes)
		# #checks convergence
		# if(CHECK_CONV(Z_ext, U_list)){
			# converged <- TRUE
		# }else{
			# curr_iter <- curr_iter + 1
			# }
	# }
	# #end of main loop
	
	# #put together return list, and returns
	# fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	# norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	# retL <- list(Z_ext=Z_ext, U=U_list, conv=converged, norm_percent = norm_percent, resids=fnorm_resids)
	# invisible(retL)	
# }

