library(STLM)
library(geoR)
library(mvtnorm)
library(ggplot2)

# This version tests posterior inference only on shorter range variograms

# Simulate a base dataset
X <- matrix(seq(0, 10, by=0.10), ncol=1)
n <- length(X)
a <- 8
b <- 2
R <- 3
eps <- 1e-4

set.seed(6)
Yobs <- drawGammaRF(X, a, b, eps, R)

#### ABC for model parameters a, b, & R ####
v     <- variog(coords=matrix(c(rep(1, n), X), nrow=n), data=Yobs, breaks=seq(0,10,by=1))$v
s_obs <- c(mean(Yobs), sd(Yobs), quantile(Yobs, c(0.05, 0.95)), v[1:3])

n_particle <- 1000
n_s <- length(s_obs)
alpha <- 0.2

prior <- matrix(10*runif(3*n_particle), ncol=3)

max_iter  <- 10
iter      <- 1
p         <- 3 # Number of model parameters

h         <- rep(Inf, max_iter) # 
d         <- matrix(0, nrow=n_particle, ncol=max_iter) # 
log_w     <- matrix(0, nrow=n_particle, ncol=max_iter)
ss_Kernel <- matrix(0, nrow=n_particle, ncol=max_iter)

posterior <- array(0, dim=c(n_particle, p, max_iter))
s_sim     <- array(0, dim=c(n_particle, n_s, max_iter))
Sigma     <- array(0, dim=c(p, p, max_iter))
Sigma_X   <- array(0, dim=c(n_s, n_s, max_iter))

omega <- matrix(1, nrow=n_s, ncol=max_iter)

start_time <- Sys.time()

# Initial iteration
for(i in 1:n_particle){
  y_sim <- drawGammaRF(X, prior[i, 1], prior[i, 2], eps, prior[i, 3])
  s_sim[i, , iter] <- c(mean(y_sim), 
                        sd(y_sim), 
                        quantile(y_sim, c(0.05, 0.95)),
                        variog(coords=matrix(c(rep(1, n), X), nrow=n), 
                               data=y_sim, 
                               breaks=seq(0,10,by=1), messages=FALSE)$v[1:3])
  d[i, iter] <- sqrt(sum((omega[, iter]*(s_obs-s_sim[i, , iter]))^2))
}
posterior[, , iter] <- prior
Sigma[, , iter]   <- diag(diag(2*cov.wt(posterior[, , iter], wt=exp(log_w[, iter]))$cov))
#Sigma_X[, , iter] <- 2*cov.wt(s_sim[, , iter], wt=exp(log_w[, iter]))$cov

log_w[, iter] <- log(1/n_particle)
h[iter + 1] <- quantile(d[, iter], alpha)
omega[, iter+1] <- 1 / apply(abs(sweep(s_sim[,, iter], 2, FUN="-", s_obs)), 2, median)

# Loop through additional iterations
eta <- rep(0, n_particle)
logPrior <- log(10^-3)

for(iter in 2:max_iter){
  # for(i in 1:n_particle){
  #   ss_Kernel[i, iter-1] <- dmvnorm(s_obs, s_sim[i, , iter-1], Sigma_X[, , iter-1], log=TRUE)
  # }
  # 
  n_accept  <- 0
  num_draws <- 0
  # abs_dev <- matrix(0, nrow=5000, ncol=n_s)
  draws <- matrix(0, nrow=10000, ncol = n_s)
  
  maxLogImpWt   <- max(log_w[, iter-1])
  # maxSSKernel   <- max(ss_Kernel[, iter-1])
  # sampleWeights <- exp(log_w[, iter-1] + ss_Kernel[, iter-1] - maxLogImpWt - maxSSKernel)
  sampleWeights <- exp(log_w[, iter-1] - maxLogImpWt)
  sampleWeights <- sampleWeights / sum(sampleWeights)
  
  while(n_accept < n_particle){
    num_draws <- num_draws + 1
    
    index <- sample(1:n_particle, 1, sampleWeights, replace=TRUE)
    # theta_star <- rmvnorm(1, mean=posterior[index, , iter-1], sigma=Sigma[, , iter-1])
    theta_star <- rmvt(1, delta=posterior[index, , iter-1], sigma=Sigma[, , iter-1], df=4)
    while(min(theta_star) < 0 | max(theta_star) > 10){
      # theta_star <- rmvnorm(1, mean=posterior[index, , iter-1], sigma=Sigma[, , iter-1])
      theta_star <- rmvt(1, delta=posterior[index, , iter-1], sigma=Sigma[, , iter-1], df=4)
    }
    
    y_sim_star <- drawGammaRF(X, theta_star[1], theta_star[2], eps, theta_star[3])
    s_sim_star <- 
      c(mean(y_sim_star), 
        sd(y_sim_star), 
        quantile(y_sim_star, c(0.05, 0.95)),
        variog(coords=matrix(c(rep(1, n), X), nrow=n), 
               data=y_sim_star, 
               breaks=seq(0,10,by=1), messages=FALSE)$v[1:3])
    
    d_star <- sqrt(sum((omega[, iter]*(s_obs-s_sim_star))^2))
    
    if(num_draws < 10000){
      # abs_dev[num_draws, ] <- abs(s_obs - s_sim_star)
      draws[num_draws, ] <- s_sim_star
    }
    
    if(d_star < h[iter]){
      n_accept <- n_accept + 1
      posterior[n_accept, , iter] <- theta_star
      s_sim[n_accept, , iter]     <- s_sim_star
      d[n_accept, iter]           <- d_star
      
      for(j in 1:n_particle){
        # Need to update importance weight for truncation of prior distribution
        eta[j] <- sampleWeights[j] + 
          # dmvnorm(theta_star, posterior[j, , iter-1], Sigma[, , iter-1], log=TRUE)
          # dmvt(theta_star, delta=posterior[j, , iter-1], sigma=Sigma[, , iter-1], log=TRUE) 
          dt((theta_star[1]-posterior[j, 1, iter-1])/sqrt(Sigma[1,1,iter-1]), df=4, log=TRUE) + 
          dt((theta_star[2]-posterior[j, 2, iter-1])/sqrt(Sigma[2,2,iter-1]), df=4, log=TRUE) + 
          dt((theta_star[3]-posterior[j, 3, iter-1])/sqrt(Sigma[3,3,iter-1]), df=4, log=TRUE) - 
          log(pt(q=(10-posterior[j, 1, iter-1])/sqrt(Sigma[1,1,iter-1]), df=4) - 
              pt(q=(0 -posterior[j, 1, iter-1])/sqrt(Sigma[1,1,iter-1]), df=4)) -
          log(pt(q=(10-posterior[j, 2, iter-1])/sqrt(Sigma[2,2,iter-1]), df=4) - 
                pt(q=(0 -posterior[j, 2, iter-1])/sqrt(Sigma[2,2,iter-1]), df=4)) -
          log(pt(q=(10-posterior[j, 3, iter-1])/sqrt(Sigma[3,3,iter-1]), df=4) - 
                pt(q=(0 -posterior[j, 3, iter-1])/sqrt(Sigma[3,3,iter-1]), df=4))
      }
      log_w[n_accept, iter] <- logPrior - max(eta) - log(sum(exp(eta-max(eta))))
    }
  }
  log_w[, iter]     <- log(exp(log_w[, iter]) / sum(exp(log_w[, iter])))
  Sigma[, , iter]   <- diag(diag(2*cov.wt(posterior[, , iter], wt=exp(log_w[, iter]))$cov))
  #Sigma_X[, , iter] <- 2*cov.wt(s_sim[, , iter], wt=exp(log_w[, iter]))$cov
  h[iter + 1] <- quantile(d[, iter], alpha)
  if(num_draws < 10000){
    omega[, iter+1] <- 1 / apply(abs(sweep(draws[1:num_draws, ], 2, FUN="-", s_obs)), 2, median)
    # omega[, iter+1] <- 1 / apply(abs_dev[1:num_draws, ], 2, median)    
  }else {
    omega[, iter+1] <- 1 / apply(abs(sweep(draws[1:10000, ], 2, FUN="-", s_obs)), 2, median)
    # omega[, iter+1] <- 1 / apply(abs_dev[1:5000, ], 2, median)
  }
  
}
end_time <- Sys.time()

# Note stopped during 8th iteration

plot_ind <- iter-1

ggplot(data=data.frame(x=X, y=Yobs), aes(x=x, y=y)) + geom_point()

ggplot(data=data.frame(a=posterior[,1,plot_ind], 
                       b=posterior[,2,plot_ind],
                       R=posterior[,3,plot_ind],
                       w=exp(log_w[,plot_ind]))) + 
  geom_point(aes(x=a,y=b, alpha=w)) + 
  theme_bw()

ggplot(data=data.frame(a=posterior[,1,plot_ind], L
                       b=posterior[,2,plot_ind],
                       R=posterior[,3,plot_ind],
                       w=exp(log_w[,plot_ind]))) + 
  geom_point(aes(x=a,y=R, alpha=w)) + 
  theme_bw()

ggplot(data=data.frame(a=posterior[,1,plot_ind], 
                       b=posterior[,2,plot_ind],
                       R=posterior[,3,plot_ind],
                       w=exp(log_w[,plot_ind]))) + 
  geom_density(aes(x=a,y=..density.., weight=w), fill="blue", alpha=0.3) + 
  theme_bw()

ggplot(data=data.frame(a=posterior[,1,plot_ind], 
                       b=posterior[,2,plot_ind],
                       R=posterior[,3,plot_ind],
                       w=exp(log_w[,plot_ind]))) + 
  geom_density(aes(x=b,y=..density.., weight=w), fill="blue", alpha=0.3) + 
  theme_bw()

ggplot(data=data.frame(a=posterior[,1,plot_ind], 
                       b=posterior[,2,plot_ind],
                       R=posterior[,3,plot_ind],
                       w=exp(log_w[,plot_ind]))) + 
  geom_density(aes(x=R,y=..density.., weight=w), fill="blue", alpha=0.3) + 
  theme_bw()
