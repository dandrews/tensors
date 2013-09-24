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

setGeneric(name="fnorm",
def=function(x){standardGeneric("fnorm")})

setGeneric(name="inner_prod",
def=function(e1,e2){standardGeneric("inner_prod")})


#### Need refClass?
# setGeneric(name="setFiber",
# def=function(x,...){standardGeneric("setFiber")})
# setGeneric(name="setSlice",
# def=function(x,...){standardGeneric("setSlice")})
# setGeneric(name="setSubtensor",
# def=function(x,...){standardGeneric("setSubtensor")})


