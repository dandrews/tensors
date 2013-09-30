#####(rTensor) Tensor Algebra and Statistical Models####

source("/Users/jamesyili/tensors/rTensor_Classes.R")
source("/Users/jamesyili/tensors/rTensor_Generics.R")
source("/Users/jamesyili/tensors/rTensor_BasicMethods.R")
source("/Users/jamesyili/tensors/rTensor_MatrixFunctions.R")
source("/Users/jamesyili/tensors/rTensor_SpecialTensors.R")
source("/Users/jamesyili/tensors/rTensor_Unfoldings.R")
source("/Users/jamesyili/tensors/rTensor_Decompositions.R")

###Main (Loading and Testing)

arr <- array(1:prod(2:7),dim=c(2,3,4,5,6,7))
arr <- array(rnorm(prod(2:7)),dim=c(2,3,4,5,6,7))
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


mat<-unfold(arrT,rs=c(4,6,2),cs=c(5,3,1))
arrT2<-fold(mat,rs=c(4,6,2),cs=c(5,3,1),modes=2:7)
identical(arrT,arrT2)

mat<-unfold(arrT,rs=c(5,6,4),cs=c(2,1,3))
arrT2<-fold(mat,rs=c(5,6,4),cs=c(2,1,3),modes=2:7)
identical(arrT,arrT2)

arr3d <- array(1:27L, dim=c(3,3,3))
arr3d <- array(rnorm(27),dim=c(3,3,3))
arr3dT <- as.tensor(arr3d)

sum1 <- 0
for(i in 1:2){
	sum1 <- sum1 + arr[i,,,,,]
}
sum2 <- 0
for(i in 1:3){
	sum2 <- sum2 + arr[,i,,,,]
}
sum3 <- 0
for(i in 1:4){
	sum3 <- sum3 + arr[,,i,,,]
}
sum4 <- 0
for(i in 1:5){
	sum4 <- sum4 + arr[,,,i,,]
}
sum5 <- 0
for(i in 1:6){
	sum5 <- sum5 + arr[,,,,i,]
}
sum6 <- 0
for(i in 1:7){
	sum6 <- sum6 + arr[,,,,,i]
}
