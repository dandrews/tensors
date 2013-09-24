#####(rTensor) Tensor Algebra and Statistical Models####

#####Various Tensor Decompositions

###(Truncated-)Higher-order SVD
#Input: Tensor tnsr, ranks = k_1 ,..., k_M (optional)
#Output: List containing a core tensor Z (all-orthogonal), and a list of matrices of sizes I_m x k_M if ranks provided, o/w I_m x I_m
hosvd <- function(tnsr,ranks=NULL){
	stopifnot(is(tnsr,"Tensor"))
	num_modes <- getNumModes(tnsr)
	
	#no truncation if ranks not provided
	if(is.null(ranks)){
		cat("!ranks not provided so left singular matrices will not be truncated.\n")
		ranks <- getModes(tnsr)	
	}

	#loops through and performs SVD on mode-m matricization of tnsr
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		U_list[[m]] <- (svd(m_unfold(tnsr,m=m))$u)[,1:ranks[m]]
	}
	
	#computes the core tensor
	Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
	
	#put together the return list, and returns
	retL <- list(Z=Z,U=U_list)
	retL
}


###CP Decomp ALS
#Input: Tensor, number of components to decompose the tensor (k)
#Output: List containing weights lambda (vector of length k), and a list of matrices U's, were each U_m is of size I_m x k
cp_als <- function(tnsr, num_components=NULL,max_iter=500, tol=1e-6){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))

	#initialization via truncated hosvd
	num_modes <- getNumModes(tnsr)
	modes <- getModes(tnsr)
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 1:num_modes){
		unfolded_mat[[m]] <- m_unfold(tnsr,m=m)
		U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
	}
	est <- tnsr
	curr_iter <- 1
	converged <- FALSE

	#set up convergence check
	fnorm_resids <- rep(0, max_iter)
	CHECK_CONV <- function(est){
		curr_resid <- fnorm(tnsr - est)
		cat("residual: ",curr_resid,"\n")
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}	

	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
		cat("iteration: ",curr_iter,"\t")
		for(m in 1:num_modes){
			V <- hamadard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
			V_inv <- solve(V)			
			tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],rev=TRUE)%*%V_inv
			lambdas <- apply(tmp,2,norm_vec)
			U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
			Z <- diagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
			est <- ttl(Z,U_list,ms=1:num_modes)
		}
		print(U_list[[1]])

		#checks convergence
		if(CHECK_CONV(est)){
			converged <- TRUE
		}else{
			curr_iter <- curr_iter + 1
			 }
	}
	#end of main loop

	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	retL <- list(lambdas=lambdas, U=U_list, conv=converged, norm_percent=norm_percent, resids=fnorm_resids)
	invisible(retL)
}


###Tucker Decomp ALS
#Input: Tensor of size I_1,...,I_M, ranks k_1,...,k_M, where M is the number of modes to Tensor
#Output: List containing a core Tensor Z and and a list of matrices U's, were each U_m is of size I_m x k_m
tucker_als <- function(tnsr,ranks=NULL,max_iter=50,tol=1e-6){
	stopifnot(is(tnsr,"Tensor"))
	if(is.null(ranks)) stop("ranks must be specified")
	
	#initialization via truncated hosvd
	num_modes <- getNumModes(tnsr)
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		U_list[[m]] <- (svd(m_unfold(tnsr,m=m))$u)[,1:ranks[m]]
	}
	curr_iter <- 1
	converged <- FALSE

	#set up convergence check
	fnorm_resids <- rep(0, max_iter)
	CHECK_CONV <- function(Z,U_list){
		est <- ttl(Z,U_list,ms=1:num_modes)
		curr_resid <- fnorm(tnsr - est)
		cat("residual: ",curr_resid,"\n")
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}
	
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	cat("iteration: ",curr_iter,"\t")
	modes_seq <- 1:num_modes
		for(m in modes_seq){
			#core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
			#truncated SVD of X
			U_list[[m]] <- (svd(m_unfold(X,m=m))$u)[,1:ranks[m]]
		}
		#compute core tensor Z
		Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

		#checks convergence
		if(CHECK_CONV(Z, U_list)){
			converged <- TRUE
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	#end of main loop
	
	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	retL <- list(Z=Z, U=U_list, conv=converged, norm_percent = norm_percent, resids=fnorm_resids)
	invisible(retL)
}


###MPCA
#Input: Tensor of size I_1,...,I_M, where m=1 is the measurement mode, ranks k_2,...,k_M, where M is the number of modes to Tensor
#Output: List containing an extended core Tensor Z_ext and and a list of matrices U's, were each U_m is of size I_m x k_m for m = 2, ..., M
mpca_als <- function(tnsr, ranks = NULL, max_iter = 500, tol=1e-6){
	if(is.null(ranks)) stop("ranks must be specified")
	stopifnot(is(tnsr,"Tensor"))

	#initialization via hosvd of M-1 modes
	num_modes <- getNumModes(tnsr)
	stopifnot(length(ranks)==(num_modes-1))
	ranks <- c(1,ranks)
	modes <- getModes(tnsr)
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 2:num_modes){
		unfolded_mat <- m_unfold(tnsr,m=m)
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
		cat("residual: ",curr_resid,"\n")
		fnorm_resids[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resids[curr_iter-1]) < tol) return(TRUE)
		FALSE
	}
	
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	cat("iteration: ",curr_iter,"\t")
	modes_seq <- 2:num_modes
		for(m in modes_seq){
			#extended core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-c(1,m)],t),ms=modes_seq[-(m-1)])
			#truncated SVD of X
			U_list[[m]] <- (svd(m_unfold(X,m=m))$u)[,1:ranks[m]]
		}
		#compute core tensor Z_ext
		Z_ext <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)
		#checks convergence
		if(CHECK_CONV(Z_ext, U_list)){
			converged <- TRUE
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	#end of main loop
	
	#put together return list, and returns
	fnorm_resids <- fnorm_resids[fnorm_resids!=0]
	norm_percent<-1-tail(fnorm_resids,1)/fnorm(tnsr)
	retL <- list(Z_ext=Z_ext, U=U_list, conv=converged, norm_percent = norm_percent, resids=fnorm_resids)
	invisible(retL)
}

###MICA - Requires W, the mixing matrix, to be specified
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

