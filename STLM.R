E1 <- function(x, alp=1e-9){return(pgamma(x, alp, lower=FALSE)/alp)}
E1inv <- function(y, alp=1e-9){return(qgamma(alp*y, alp, lower=FALSE))}

cyl.corr <- function(d, R){
  return((2/pi)*acos(d/(2*R)) - (d/(2*pi*R^2))*sqrt(4*R^2-d^2))
}

plot.circle <- function(c, r, Color="blue"){
  x.axis <- seq(c[1]-r, c[1]+r, by=r/100)
  lines(x.axis, c[2] + sqrt(r^2-(x.axis-c[1])^2), col=Color)
  lines(x.axis, c[2] - sqrt(r^2-(x.axis-c[1])^2), col=Color)
}

drawGammaLM <- function(eps, a, b, nc){
  cdf <- seq(1e-5, 1-1e-5, length=2500)
  x <- E1inv(E1(b*eps)*(1-cdf))
  J <- rpois(1, nc)
  output <- x[findInterval(runif(J, eps, 1), cdf)]
  return(output)
}

drawGammaRF <- function(X, a, b, eps, R, sd=NULL){
  set.seed(sd)
  
  norm.constant <- a*E1(b*eps)
  
  # Storage for Levy mass "jumps"
  mass.pts <- matrix(0, nrow=10*round(norm.constant*nrow(X)), ncol=3)
  num.added <- rep(0, nrow(X))
  total.mass.pts <- 0
  
  # Draw mass points for first location
  levy.draws <- drawGammaLM(1e-5, 10, 1, norm.constant)
  num.draws <- length(levy.draws)
  draw.radius <- sqrt(runif(num.draws, 0, R))
  draw.angle <- runif(num.draws, 0, 2*pi)
  
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 1] <- levy.draws
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 2] <- draw.radius*cos(draw.angle) + X[1, 1]
  mass.pts[(total.mass.pts+1):(total.mass.pts+num.draws), 3] <- draw.radius*sin(draw.angle) + X[1, 2]
  
  num.added[1] <- num.draws
  total.mass.pts <- total.mass.pts + num.draws
  
  for(i in 2:nrow(X)){
    # Draw mass points for other locations, but discard points that fall 
    # in the regions of previous locations
    levy.draws <- drawGammaLM(1e-5, 10, 1, norm.constant)
    num.draws <- length(levy.draws)
    draw.radius <- sqrt(runif(num.draws, 0, R))
    draw.angle <- runif(num.draws, 0, 2*pi)
    
    x.loc <- X[i,1] + draw.radius*cos(draw.angle)
    y.loc <- X[i,2] + draw.radius*sin(draw.angle)
    
    valid.new.pt <- rep(TRUE, num.draws)
    
    for(j in 1:(i-1)){
      valid.new.pt <- valid.new.pt & !(((x.loc-X[j,1])^2 + (y.loc-X[j,2])^2) < R^2)
    }
    
    num.new <- sum(valid.new.pt)
    num.added[i] <- num.new
    
    if(num.new > 0){
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[valid.new.pt]
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- x.loc[valid.new.pt]
      mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- y.loc[valid.new.pt]
      total.mass.pts <- total.mass.pts + num.new
    }
  }
  Y <- rep(0, nrow(X))
  for(i in 1:nrow(X)){
    Y[i] <- sum(mass.pts[(mass.pts[,2]-X[i,1])^2 + (mass.pts[,3]-X[i,2])^2 < R^2,1])
  }
  
  return(cbind(X, Y, num.added))
}

# Marginal gamma distribution parameters
a <- 10
b <- 1
eps <- 1e-5

# Cylinder based radius
R <- 0.3

# Grid on which to evaluate data
X <- expand.grid(seq(0.05, 0.95, length=5), seq(0.05, 0.95, length=5))

# test <- drawGammaRF(10, 1, 1e-5, R, 0.05, 0.95, 0.05, 0.95, 20, 1)

marg.results <- matrix(0, nrow=1000, ncol=25)
for(iter in 1:1000){
  junk <- drawGammaRF(X, 10, 1, 1e-5, R=R, sd=iter+1)
  marg.results[iter,] <- junk$Y
}

ggplot(data=data.frame(expand.grid(seq(0.05, 0.95, length=5), seq(0.05, 0.95, length=5)), Y=marg.results[1000,]), 
       aes(x=Var1, y=Var2)) + 
  geom_raster(aes(fill=Y))

X <- expand.grid(seq(0.05, 0.95, length=20), seq(0.05, 0.95, length=20))

plot(x.loc, y.loc, xlim=c(-0.3,1.3), ylim=c(-0.3, 1.3))
plot.circle(unlist(X[1,]), 0.3)
plot.circle(unlist(X[2,]), 0.3)

plot(mass.pts[1:102,2:3], col="green", pch=16, xlim=c(-0.3,1.3), ylim=c(-0.3, 1.3))
lines(mass.pts[103:138,2:3], col="red", pch=16, type="p")
plot.circle(unlist(X[1,]), 0.3)
plot.circle(unlist(X[2,]), 0.3)


hist(marg.results[,3], freq=FALSE)
lines(seq(0.01, 25, by=0.01), dgamma(seq(0.01, 25, by=0.01), 10, 1), col="red")

cor(marg.results[,4], marg.results[,5])
cyl.corr(0.225, 0.3)


X <- expand.grid(seq(0.25, 0.75, length=2), seq(0.25, 0.75, length=2))
drawGammaRF(X, 10, 1, 1e-5, R=0.5, sd=6)


ggplot(data=data.frame(expand.grid(seq(0.05, 0.95, length=5), seq(0.05, 0.95, length=5)), Y=drawGammaRF(X, 10, 1, 1e-5, R=0.001, sd=6)$Y), 
       aes(x=Var1, y=Var2)) + 
  geom_raster(aes(fill=Y))
