#####(rTensor) Tensor Algebra and Statistical Models####
#####Special Tensors
###Create a Random Tensor
rand_tensor <- function(modes=c(3,4,5)){
	as.tensor(array(runif(prod(modes)), dim=modes))
}
###Create a Superdiagonal Tensor
superdiagonal_tensor <- function(num_modes,len,elements=1L){
	modes <- rep(len,num_modes)
	arr <- array(0, dim = modes)
	if(length(elements)==1) elements <- rep(elements,len)
	for (i in 1:len){
		txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
		eval(parse(text=txt))
	}
	as.tensor(arr)
}
###Create an 3d Identity Tensor
identity_tensor3d <- function(modes){
	if(length(modes)!=3L) stop("identity tensor only implemented for 3d so far")
	n <- modes[1]
	stopifnot(n==modes[2])
	arr <- array(0,dim=modes)
	arr[,,1] <- diag(1,n,n)
	as.tensor(arr)
}
