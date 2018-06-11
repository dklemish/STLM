# ABC.1D <- function(X, y, prior, alpha, nparticle, nbins, 
#                    corrStruct="Cylinder", sd=NULL,
#                    maxDist, eps=1e-4, maxIter){
#   set.seed(sd)
#   
#   if(corrStruct=="Cylinder"){
#     type <- 1
#   }else if(corrStruct=="Cone"){
#     type <- 2
#   }
#   
#   distMatrix <- as.matrix(dist(X))
#   
#   prior.param <- unlist(prior)
#   
#   #result <- .Call('_STLM_ABC_1D', PACKAGE='STLM', X, shape, rate, eps, R, type, distMatrix)
#   
#   return(prior.param)
# }