#E1 <- function(x, alp=1e-9){return(pgamma(x, alp, lower=FALSE)/alp)}
#E1inv <- function(y, alp=1e-9){return(qgamma(alp*y, alp, lower=FALSE))}

cyl.corr <- Vectorize(function(d, R, dim){
  if(d > 2*R){
    return(0)
  }else if(dim == 1){
    return(1-0.5*d/R)
  }else if(dim == 2){
    return((2*R^2*acos(d/(2*R))-(0.5*d*sqrt(4*R^2 - d^2))) / (pi*R^2))
  }
})

cone.corr <- Vectorize(function(d, R, dim){
  if(d > 2*R){
    return(0)
  }else if(dim == 1){
    return(0.25*(2-d/R)^2)
  }else if(dim == 2){
    #return((2/pi)*acos(d/(2*R)) - (d/(2*pi*R^2))*sqrt(4*R^2-d^2))
    return("Not correct yet!")
  }
})

gauss.corr <- Vectorize(function(d, R, sig, dim){
  if(d > 2*R){
    return(0)
  }else if(dim == 1){
    return((pnorm(R, 0, sig2) - pnorm(d/2, 0, sig2)) / (pnorm(R, 0, sig2) - 0.5))
  }else if(dim == 2){
    #return((2/pi)*acos(d/(2*R)) - (d/(2*pi*R^2))*sqrt(4*R^2-d^2))
    return("Not correct yet!")
  }
})

drawGammaRF <- function(X, shape, rate, eps, R=1, corrStruct="Cylinder", sd=NULL){
  set.seed(sd)
  if(corrStruct=="Cylinder"){
    type <- 1
  }else if(corrStruct=="Cone"){
    type <- 2
  }
  
  distMatrix <- as.matrix(dist(X))
  
  #print(distMatrix)
  
  if(ncol(X) == 1){
    result <- .Call('_STLM_drawGammaRF_1D', PACKAGE='STLM', X, shape, rate, eps, R, type, distMatrix)  
  }else{
    print("Not implemented yet!")
    #result <- .Call('_STLM_drawGammaRF_2D', PACKAGE='STLM', X, shape, rate, eps, R, type, distMatrix)
  }
  
  return(result)
}

# drawGammaRF <- function(X, a, b, eps, R, shape="Cylinder", mult=10, sd=NULL){
#   set.seed(sd)
#   
#   norm.constant <- a*E1(b*eps)
#   
#   if(shape=="Cylinder"){
#     nSpatial <- 3 # mass value + (x, y)
#   }else if(shape=="Cone"){
#     nSpatial <- 4 # mass vale + (x, y, h)
#   }
#   
#   # Storage for Levy mass "jumps"
#   mass.pts <- matrix(0, nrow=mult*round(norm.constant*nrow(X)), ncol=nSpatial)
#   num.added <- rep(0, nrow(X))
#   total.mass.pts <- 0
#   
#   # Draw mass points for first location
#   levy.draws <- drawGammaLM(eps, a, b, norm.constant)
#   num.draws <- length(levy.draws)
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 1] <- levy.draws
#   
#   if(shape=="Cylinder"){
#     draw.radius <- sqrt(runif(num.draws, 0, R))
#     draw.angle <- runif(num.draws, 0, 2*pi)
#     mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 2] <- draw.radius*cos(draw.angle) + X[1, 1]
#     mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 3] <- draw.radius*sin(draw.angle) + X[1, 2]
#   }
#   
#   num.added[1] <- num.draws
#   total.mass.pts <- total.mass.pts + num.draws
#   
#   for(i in 2:nrow(X)){
#     # Draw mass points for other locations, but discard points that fall 
#     # in the regions of previous locations
#     levy.draws <- drawGammaLM(eps, a, b, norm.constant)
#     num.draws  <- length(levy.draws)
# 
#     if(shape=="Cylinder"){
#       draw.radius <- sqrt(runif(num.draws, 0, R))
#       draw.angle  <- runif(num.draws, 0, 2*pi)
# 
#       x.loc <- X[i,1] + draw.radius*cos(draw.angle)
#       y.loc <- X[i,2] + draw.radius*sin(draw.angle)
# 
#       valid.new.pt <- rep(TRUE, num.draws)
#       
#       # Check previous locations whether new mass positions fall in those location's shapes
#       for(j in 1:(i-1)){
#         #valid.new.pt <- valid.new.pt & !(((x.loc-X[j,1])^2 + (y.loc-X[j,2])^2) < R^2)
#         valid.new.pt <- !(((x.loc-X[j,1])^2 + (y.loc-X[j,2])^2) < R^2)
#       }
#       
#       num.new <- sum(valid.new.pt)
#       num.added[i] <- num.new
# 
#       if(num.new > 0){
#         mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[valid.new.pt]
#         mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- x.loc[valid.new.pt]
#         mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- y.loc[valid.new.pt]
#         total.mass.pts <- total.mass.pts + num.new
#       }
#     }
#   }
#   
#   # Calculate value of response at each location
#   Y <- rep(0, nrow(X))
#   
#   if(shape=="Cylinder"){
#     for(i in 1:nrow(X)){
#       Y[i] <- sum(mass.pts[(mass.pts[,2]-X[i,1])^2 + (mass.pts[,3]-X[i,2])^2 < R^2,1])
#     }    
#   }
# 
#   return(cbind(X, Y, num.added))
# }

drawGammaRF_R <- function(X, a, b, eps, R, sd=NULL){
  set.seed(sd)
  
  nSpatial <- 3 # mass value, x, h
  nLoc <- nrow(X)
  mult <- 10000
  
  # Storage for Levy mass "jumps"
  mass.pts <- matrix(0, nrow=mult*nLoc, ncol=nSpatial)
  num.added <- rep(0, nLoc)
  total.mass.pts <- 0
  
  # Draw mass points for first location
  levy.draws <- drawGamma(10, a, b, eps)
  num.draws <- sum(levy.draws > 0)
  
  ### Note: Not correct uniform distribution on triangle!
  draw.loc <- X[1] + runif(num.draws, -R, R)
  draw.height <- runif(num.draws, 0, (R - abs(draw.loc - X[1])) / R^2)
  
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 1] <- levy.draws[1:num.draws]
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 2] <- draw.loc
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 3] <- draw.height
  
  # plot(-9:1, (0:10)*0.01, type="l", xlim=c(-9,11))
  # lines(1:11, (0:10)*-0.01 + 0.1)
  # lines(draw.loc, draw.height, col="red", type="p")
  
  num.added[1] <- num.draws
  total.mass.pts <- total.mass.pts + num.draws
  
  for(i in 2:nLoc){
    # Draw mass points for other locations, but discard points that fall
    # in the regions of previous locations
    levy.draws <- drawGamma(10, a, b, eps)
    num.draws  <- sum(levy.draws > 0)
    
    draw.loc    <- X[i,1] + runif(num.draws, -R, R)
    draw.height <- runif(num.draws, 0, (R - abs(draw.loc - X[i,1])) / R^2) 
    
    valid.new.pt <- rep(TRUE, num.draws)
    
    # Check previous locations whether new mass positions fall in those location's shapes
    for(j in 1:(i-1)){
      valid.new.pt <- valid.new.pt & 
        (abs(X[j,1] - draw.loc) > R) & 
        (draw.height > (R-abs(draw.loc - X[j,1]))/R^2)
    }
    
    num.new <- sum(valid.new.pt)
    num.added[i] <- num.new
    
    if(num.new > 0){
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[1:num.draws][valid.new.pt]
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- draw.loc[valid.new.pt]
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- draw.height[valid.new.pt]
      total.mass.pts <- total.mass.pts + num.new
    }
  }
  mass.pts <- mass.pts[1:total.mass.pts, ]
  
  # Calculate value of response at each location
  Y <- rep(0, nLoc)
  
  for(i in 1:nLoc){
    Y[i] <- sum(mass.pts[(abs(mass.pts[,2]-X[i,1]) < R) & 
                           (mass.pts[,3] < (R-abs(mass.pts[,2] - X[i,1]))/R^2),
                         1])
  }
  
  
  return(Y)
}
