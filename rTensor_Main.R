#####(rTensor) Tensor Algebra and Statistical Models####

source("/Users/jamesyili/tensors/rTensor_Classes.R")
source("/Users/jamesyili/tensors/rTensor_Generics.R")
source("/Users/jamesyili/tensors/rTensor_Methods.R")
source("/Users/jamesyili/tensors/rTensor_Functions.R")
source("/Users/jamesyili/tensors/rTensor_Special.R")
source("/Users/jamesyili/tensors/rTensor_Decompositions.R")

###Main (Loading and Testing)

arr <- array(1:prod(2:7),dim=c(2,3,4,5,6,7))
mat <- matrix(1:100,nrow=10)
vec <- 1:10
arrT <- as.tensor(arr)
matT <- as.tensor(mat)
vecT <- as.tensor(vec)

testlist <- list(arrT,matT,vecT)
lapply(testlist,print)
lapply(testlist,show)
lapply(testlist,getNumModes)
lapply(testlist,getModes)
lapply(testlist,getModenames)
lapply(testlist,getData)

