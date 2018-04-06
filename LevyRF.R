a <- 1
b <- 1

nu <- 
  Vectorize(function(x, a, b){
    return(a * exp(-b*x) / x)
  })

rgam <- function(n, a, b){
  output <- rep(0, n)
  
  num.partition <- 10
  x <- seq(0.001, qgamma(1-1e-8, a, b), length=101)  

  for(k in 1:n){
    
    GAMMA <- rep(0, num.partition)
    integral <- rep(0, num.partition)
    
    for(i in 1:num.partition){
      v_j <- integrate(nu, x[i], x[i+1], a, b)
      n_j <- rpois(1, v_j$value)
      
      s <- seq(x[i], x[i+1], length=101)
      integral[1] <- (nu(s[2], a, b) + nu(s[1], a, b)) / (s[2] - s[1])
      for(j in 2:100){
        integral[j] <- integral[j-1] + (nu(s[j+1], a, b) + nu(s[j], a, b)) / (s[j+1] - s[j])
      }
      integral <- integral / integral[100]
      
      U <- runif(n_j, 0, 1)
      
      if(n_j > 0){
        for(j in 1:n_j){
          GAMMA[i] <- GAMMA[i] + s[which.min(U[j] > integral)]
        }
      }
    }
    output[k] <- sum(GAMMA)
  }
  return(output)
}
