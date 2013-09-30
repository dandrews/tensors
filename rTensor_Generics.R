#####(rTensor) Tensor Algebra and Statistical Models####

#####Generic Definitions
setGeneric(name="getModes",
def=function(x){standardGeneric("getModes")})

setGeneric(name="getNumModes",
def=function(x){standardGeneric("getNumModes")})

setGeneric(name="getModenames",
def=function(x){standardGeneric("getModenames")})

setGeneric(name="getData",
def=function(x){standardGeneric("getData")})

setGeneric(name="getFiber",
def=function(x,...){standardGeneric("getFiber")})

setGeneric(name="getSlice",
def=function(x,...){standardGeneric("getSlice")})

setGeneric(name="getSubtensor",
def=function(x,...){standardGeneric("getSubtensor")})

setGeneric(name="rs_unfold",
def=function(x,...){standardGeneric("rs_unfold")})

setGeneric(name="cs_unfold",
def=function(x,...){standardGeneric("cs_unfold")})

setGeneric(name="unfold",
def=function(x,...){standardGeneric("unfold")})

setGeneric(name="rs_fold",
def=function(x,...){standardGeneric("rs_fold")})

setGeneric(name="cs_fold",
def=function(x,...){standardGeneric("cs_fold")})

setGeneric(name="fold",
def=function(x,...){standardGeneric("fold")})

setGeneric(name="modeSum",
def=function(x,...){standardGeneric("modeSum")})

setGeneric(name="fnorm",
def=function(x){standardGeneric("fnorm")})

setGeneric(name="inner_prod",
def=function(x1,x2){standardGeneric("inner_prod")})

setGeneric(name="hosvd",
def=function(x,...){standardGeneric("hosvd")})

setGeneric(name="cp_als",
def=function(x,...){standardGeneric("cp_als")})

setGeneric(name="tucker_als",
def=function(x,...){standardGeneric("tucker_als")})

setGeneric(name="mpca_als",
def=function(x,...){standardGeneric("mpca_als")})

setGeneric(name="fft_svd",
def=function(x,...){standardGeneric("fft_svd")})














#### Need refClass?
# setGeneric(name="setFiber",
# def=function(x,...){standardGeneric("setFiber")})
# setGeneric(name="setSlice",
# def=function(x,...){standardGeneric("setSlice")})
# setGeneric(name="setSubtensor",
# def=function(x,...){standardGeneric("setSubtensor")})


