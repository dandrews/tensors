#####(rTensor) Tensor Algebra and Statistical Models####

#####Special Tensors

###Diagonal Tensor
diagonal_tensor <- function(num_modes,len,elements=1L){
	modes <- rep(len,num_modes)
	arr <- array(0, dim = modes)
	if(length(elements)==1) elements <- rep(elements,len)
	for (i in 1:len){
		txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
		eval(parse(text=txt))
	}
	tnsr <- new("ndTensor", data=arr)
	tnsr
}



