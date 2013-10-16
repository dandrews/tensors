#Timing and Percent norm recovered
source("/Users/jamesyili/tensors/rTensor_Classes.R")
source("/Users/jamesyili/tensors/rTensor_Generics.R")
source("/Users/jamesyili/tensors/rTensor_BasicMethods.R")
source("/Users/jamesyili/tensors/rTensor_MatrixFunctions.R")
source("/Users/jamesyili/tensors/rTensor_SpecialTensors.R")
source("/Users/jamesyili/tensors/rTensor_Unfoldings.R")
source("/Users/jamesyili/tensors/rTensor_Decompositions.R")

#tic toc functions
tic <- function (gcFirst = TRUE,overwrite=TRUE) {
   if(gcFirst) gc(FALSE)
   tic <- proc.time()
   ticExists <- ".tic"%in%ls(all.names=TRUE,envir=baseenv())
   if(overwrite||!ticExists){
   	assign(".tic", tic, envir=baseenv())
   	}
   	else{
   		stop("Another timing function running")
   		}
   invisible(tic)
}

toc <- function (kill=TRUE,pr=FALSE) {
   toc <- proc.time()
   tic <- get(".tic", envir=baseenv())
   if(pr) print(toc - tic)
   if(kill) rm(.tic, envir=baseenv())
   invisible(toc - tic)
}

#1000 x 1000 x 100
tnsr <- new("Tensor",3L,c(1000L,1000L,100L),data=runif(1e08))
object.size(tnsr) #800 001 296 bytes (800 MB)

save.image(file="firstcomparisons.RData")

#hosvd
tic()
hosvdD <- hosvd(tnsr)
time.hosvdD <- toc()

save.image(file="firstcomparisons.RData")

#cp
tic()
cpD <- cp_als(tnsr)
time.cpD <- toc()

save.image(file="firstcomparisons.RData")

#tucker
tic()
tuckerD <- tucker_als(tnsr)
time.tukcerD <- toc()

save.image(file="firstcomparisons.RData")

#mpca
tic()
mpcaD <- mpca_als(tnsr)
time.mpcaD <- toc()

save.image(file="firstcomparisons.RData")

#t_svd
tic()
tsvdD <- t_svd(tnsr)
time.tsvdD <- toc()

save.image(file="firstcomparisons.RData")
