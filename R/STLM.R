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
    return((1-d/(2*R))^2)
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
  
  # Calculate distances of all observed points
  distMatrix <- as.matrix(dist(X))
  
  if(ncol(X) == 1){
    result <- .Call('_STLM_drawGammaRF_1D', PACKAGE='STLM', X, shape, rate, eps, R, type, distMatrix)  
  }else{
    result <- .Call('_STLM_drawGammaRF_2D', PACKAGE='STLM', X, shape, rate, eps, R, type, distMatrix)
  }
  
  return(result)
}

cone_height1D <- function(d, R){
    return((R-abs(d))/R^2)
}

# drawGammaRF_1D_Cone_R <- function(X, a, b, eps, R, sd=NULL){
#   set.seed(sd)
#   
#   nSpatial <- 3 # mass value, x, h
#   nLoc <- nrow(X)
#   mult <- 1000
#   
#   distMatrix <- as.matrix(dist(X))
#   
#   # Storage for Levy mass "jumps"
#   mass.pts <- matrix(0, nrow=mult*nLoc, ncol=nSpatial)
#   num.added <- rep(0, nLoc)
#   total.mass.pts <- 0
#   
#   # Draw mass points for first location
#   levy.draws <- drawGamma(10, a, b, eps)
#   num.draws <- sum(levy.draws > 0)
#   num.draws.pos <- rbinom(num.draws, 1, 0.5)
#   
#   draw.loc    <- rep(0, num.draws)
#   draw.loc[num.draws.pos == 0] <- (X[1] - R) + R*rbeta(num.draws - sum(num.draws.pos), 2, 1)
#   draw.loc[num.draws.pos == 1] <- X[1] + R*rbeta(sum(num.draws.pos), 1, 2)
#   draw.height <- runif(num.draws, 0, cone_height1D(draw.loc - X[1], R))
#   
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 1] <- levy.draws[1:num.draws]
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 2] <- draw.loc
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 3] <- draw.height
#   
#   num.added[1] <- num.draws
#   total.mass.pts <- total.mass.pts + num.draws
#   
#   for(i in 2:nLoc){
#     # Draw mass points for other locations, but discard points that fall
#     # in the regions of previous locations
#     levy.draws <- drawGamma(10, a, b, eps)
#     num.draws  <- sum(levy.draws > 0)
#     draw.loc    <- rep(0, num.draws)
#     
#     num.draws.pos <- rbinom(num.draws, 1, 0.5)
#     
#     draw.loc[num.draws.pos == 0] <- (X[i] - R) + R*rbeta(num.draws - sum(num.draws.pos), 2, 1)
#     draw.loc[num.draws.pos == 1] <- X[i] + R*rbeta(sum(num.draws.pos), 1, 2)
#     draw.height <- runif(num.draws, 0, cone_height1D(draw.loc - X[i], R)) 
#     
#     valid.new.pt <- rep(TRUE, num.draws)
#     
#     # Check previous locations whether new mass positions fall in those location's shapes
#     for(j in 1:(i-1)){
#       if(distMatrix[i,j] < 2*R){
#         valid.new.pt <- valid.new.pt & 
#           (draw.height > cone_height1D(draw.loc - X[j], R))
#       }
#     }
#     
#     num.new <- sum(valid.new.pt)
#     num.added[i] <- num.new
#     
#     if(num.new > 0){
#       mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[1:num.draws][valid.new.pt]
#       mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- draw.loc[valid.new.pt]
#       mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- draw.height[valid.new.pt]
#       total.mass.pts <- total.mass.pts + num.new
#     }
#   }
#   mass.pts <- mass.pts[1:total.mass.pts, ]
#   
#   # Calculate value of response at each location
#   Y <- rep(0, nLoc)
#   
#   for(i in 1:nLoc){
#     Y[i] <- sum(mass.pts[mass.pts[,3] < cone_height1D(X[i] - mass.pts[,2], R), 1])
#   }
#   
#   return(Y)
# }
# 
# drawGammaRF_2D_Cyl_R <- function(X, a, b, eps, R, sd=NULL){
#   set.seed(sd)
#   
#   distMatrix <- as.matrix(dist(X))  
#   
#   nSpatial <- 4 # mass value, x, y, h
#   nLoc <- nrow(X)
#   mult <- 1000
#   
#   # Storage for Levy mass "jumps"
#   mass.pts <- matrix(0, nrow=mult*nLoc, ncol=nSpatial)
#   num.added <- rep(0, nLoc)
#   total.mass.pts <- 0
#   
#   # Draw mass points for first location
#   levy.draws <- drawGamma(10, a, b, eps)
#   num.draws <- sum(levy.draws > 0)
#   
#   draw.radius <- sqrt(runif(num.draws, 0, R))
#   draw.angle <- runif(num.draws, 0, 2*pi)
#   draw.loc <- matrix(0, nrow=num.draws, ncol=2)
#   draw.loc[,1] <- X[1,1] + draw.radius*cos(draw.angle)
#   draw.loc[,2] <- X[1,2] + draw.radius*sin(draw.angle)
#   
#   #draw.height <- runif(num.draws, 0, (R - abs(draw.loc - X[1])) / R^2)
#   
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 1] <- levy.draws[1:num.draws]
#   mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 2:3] <- draw.loc
#   #mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 3] <- draw.height
#   
#   num.added[1] <- num.draws
#   total.mass.pts <- total.mass.pts + num.draws
#   
#   for(i in 2:nLoc){
#     # Draw mass points for other locations, but discard points that fall
#     # in the regions of previous locations
#     levy.draws <- drawGamma(10, a, b, eps)
#     num.draws  <- sum(levy.draws > 0)
#     draw.radius <- sqrt(runif(num.draws, 0, R))
#     draw.angle <- runif(num.draws, 0, 2*pi)
#     draw.loc <- matrix(0, nrow=num.draws, ncol=2)
#     draw.loc[,1] <- X[i,1] + draw.radius*cos(draw.angle)
#     draw.loc[,2] <- X[i,2] + draw.radius*sin(draw.angle)
#     #draw.height <- runif(num.draws, 0, (R - abs(draw.loc - X[i,1])) / R^2) 
#     
#     valid.new.pt <- rep(TRUE, num.draws)
#     
#     # Check previous locations whether new mass positions fall in those location's shapes
#     for(j in 1:(i-1)){
#       if(distMatrix[i,j] < 2*R){
#         valid.new.pt <- valid.new.pt & !(((draw.loc[,1]-X[j,1])^2 + (draw.loc[,1]-X[j,2])^2) < R^2)        
#       }
#     }
# 
#     num.new <- sum(valid.new.pt)
#     num.added[i] <- num.new
#     
#     if(num.new > 0){
#       mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[1:num.draws][valid.new.pt]
#       mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2:3] <- draw.loc[valid.new.pt]
#       #mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- draw.height[valid.new.pt]
#       total.mass.pts <- total.mass.pts + num.new
#     }
#   }
#   mass.pts <- mass.pts[1:total.mass.pts, ]
#   
#   # Calculate value of response at each location
#   Y <- rep(0, nLoc)
#   
#   for(i in 1:nLoc){
#     Y[i] <- sum(mass.pts[(abs(mass.pts[,2]-X[i,1]) < R) & 
#                            (mass.pts[,3] < (R-abs(mass.pts[,2] - X[i,1]))/R^2),
#                          1])
#   }
#   
#   return(Y)
# }
