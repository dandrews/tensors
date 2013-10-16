#####(rTensor) Tensor Algebra and Statistical Models####

source("/Users/jamesyili/tensors/rTensor_Classes.R")
source("/Users/jamesyili/tensors/rTensor_Generics.R")
source("/Users/jamesyili/tensors/rTensor_BasicMethods.R")
source("/Users/jamesyili/tensors/rTensor_MatrixFunctions.R")
source("/Users/jamesyili/tensors/rTensor_SpecialTensors.R")
source("/Users/jamesyili/tensors/rTensor_Unfoldings.R")
source("/Users/jamesyili/tensors/rTensor_Decompositions.R")

#####Compilation of Tests
###Classes
##Tensor creation
#using new
tnsr <- new("Tensor",3L,c(10L,20L,30L),letters[1:3],data=runif(6000))
#from vectors
vec <- runif(100)
vecT <- as.tensor(vec)
object.size(vec); object.size(vecT) #Really noticeable for small vectors
#from matrices
mat <- matrix(runif(1000),nrow=100,ncol=10)
matT <- as.tensor(mat)
object.size(mat); object.size(matT)
#from arrays
indices <- c(10,30,100,300)
arr <- array(runif(prod(indices)), dim = indices)
arrT <- as.tensor(arr)
object.size(arr); object.size(arrT)

###BasicMethods
tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#getters, show, print, head, tail
getModes(tnsr)
getNumModes(tnsr)
getModenames(tnsr)
getData(tnsr)
tnsr
print(tnsr)
head(tnsr)
tail(tnsr)
#element-wise operation
tnsr2 <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
tnsrsum <- tnsr + tnsr2
tnsrdiff <- tnsr - tnsr2
tnsrelemprod <- tnsr * tnsr2
tnsrelemquot <- tnsr / tnsr2
for (i in 1:3L){
	for (j in 1:4L){
		for (k in 1:5L){
			stopifnot(tnsrsum[i,j,k]==tnsr[i,j,k]+tnsr2[i,j,k])
			stopifnot(tnsrdiff[i,j,k]==tnsr[i,j,k]-tnsr2[i,j,k])
			stopifnot(tnsrelemprod[i,j,k]==tnsr[i,j,k]*tnsr2[i,j,k])
			stopifnot(tnsrelemquot[i,j,k]==tnsr[i,j,k]/tnsr2[i,j,k])
		}
	}
}
#subsetting
tnsr[1,2,3]
tnsr[3,1,]
tnsr[,,5]
#modeSum
modeSum(tnsr,3)
modeMean(tnsr,1)
#sweep
sweep(tnsr,m=c(2,3),stat=1,func='-')
sweep(tnsr,m=1,stat=10,func='/')
#fnorm
fnorm(tnsr)
#inner product
innerProd(tnsr,tnsr2)

###Unfoldings
#unfolds
matT1<-cs_unfold(tnsr,m=3)
matT2<-rs_unfold(tnsr,m=2)
identical(matT1,unfold(tnsr,rs=c(1,2),cs=c(3)))
identical(matT2,unfold(tnsr,rs=2,cs=c(1,3)))
matT3<-unfold(tnsr,rs=2,cs=c(3,1))
#folds
identical(cs_fold(matT1,m=3,modes=c(3,4,5)),tnsr)
identical(rs_fold(matT2,m=2,modes=c(3,4,5)),tnsr)
identical(fold(matT3,rs=2,cs=c(3,1),modes=c(3,4,5)),tnsr)

###Operations
tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#ttm
mat <- matrix(runif(50),ncol=5)
ttm(tnsr,mat,m=3)
#ttv
vec <- runif(4)
ttv(tnsr,vec,m=2)
#ttl
lizt <- list('mat1' = matrix(runif(30),ncol=3), 'mat2' = matrix(runif(40),ncol=4),'mat3' = matrix(runif(50),ncol=5))
ttl(tnsr,lizt,ms=c(1,2,3))
#t
identical(t(tnsr)@data[,,1],t(tnsr@data[,,1]))
identical(t(tnsr)@data[,,2],t(tnsr@data[,,5]))
identical(t(t(tnsr)),tnsr)
#%*%
tnsr2 <- new("Tensor",3L,c(4L,3L,5L),data=runif(60))
tnsr%*%tnsr2

###MatrixFunctions
#hamadard_list
lizt <- list('mat1' = matrix(runif(40),ncol=4), 'mat2' = matrix(runif(40),ncol=4),'mat3' = matrix(runif(40),ncol=4))
dim(hamadard_list(lizt))
#kronecker_list
smalllizt <- list('mat1' = matrix(runif(12),ncol=4), 'mat2' = matrix(runif(12),ncol=4),'mat3' = matrix(runif(12),ncol=4))
dim(kronecker_list(smalllizt))
#khartri_rao
dim(khatri_rao(matrix(runif(12),ncol=4),matrix(runif(12),ncol=4)))
#khartri_rao_list
dim(khatri_rao_list(smalllizt))
#circ_mat
circ_mat(1:10L)

###Decompositions
tnsr <- new("Tensor",3L,c(60L,70L,80L),data=runif(336000))
smalltnsr <- new("Tensor",3L,c(10L,10L,10L),data=runif(1000))
#hosvd
hosvdD <-hosvd(tnsr)
hosvdD$resid
hosvdD2 <-hosvd(tnsr,ranks=c(6L,7L,8L))
hosvdD2$resid
#cp_als
cpD <- cp_als(tnsr,num_components=30) #(30^3)/(60*70*80) = 0.08035714
cpD$conv #did not converge with 500 iterations
cpD$norm_percent # 51%
plot(cpD$resids) 
smallcpD <- cp_als(smalltnsr,num_components=5)
smallcpD$conv 
smallcpD$norm_percent # 57%
plot(smallcpD$resids)
#tucker_als
tuckerD <- tucker_als(tnsr,ranks=c(30,35,40))
tuckerD$conv #did not converge with 500 iterations
tuckerD$norm_percent #56%
plot(tuckerD$resids/fnorm(tnsr))
smalltuckerD <- tucker_als(smalltnsr,ranks=c(5,6,7))
smalltuckerD$conv
smalltuckerD$norm_percent  #63%
plot(smalltuckerD$resids)
#mpca_als
mpcaD <- mpca_als(tnsr,ranks=c(30,30))
mpcaD$conv #converged
mpcaD$norm_percent #56%
plot(mpcaD$resids)
smallmpcaD <- mpca_als(smalltnsr,ranks=c(5,5))
smallmpcaD$conv
smallmpcaD$norm_percent #63%
plot(smallmpcaD$resids)
#tsvd3d
tsvdD <- t_svd3d(tnsr)
1 - fnorm(t_svd_reconstruct(tsvdD)-tnsr)/fnorm(tnsr) #98.5%
smalltsvdD <- t_svd3d(smalltnsr)
1 - fnorm(t_svd_reconstruct(smalltsvdD)-smalltnsr)/fnorm(smalltnsr)

###SpecialTensors
#random tensor
rand_tensor()
rand_tensor(c(8,2,100,4))
#superdiagonal tensor
superdiagonal_tensor(3,4)@data
#identity tensor
identity_tensor3d(c(3,3,10))@data

