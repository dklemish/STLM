# dmvt <- function(x, mu = rep(0, nrow(Sigma)), Sigma=diag(2), df){
#   return(.Call('_STLM_dmvt', PACKAGE='STLM', mu, Sigma, df))
# }
# 
# dmvnorm <- function(x, mu = rep(0, nrow(Sigma)), Sigma=diag(2)){
#   return(.Call('_STLM_dmvnorm', PACKAGE='STLM', mu, Sigma))
# }
# 
# dmvt.chol <- function(x, mu = rep(0, nrow(L)), L=diag(2), df){
#   return(.Call('_STLM_dmvt_chol', PACKAGE='STLM', mu, L, df))
# }
# 
# dmvnorm.chol <- function(x, mu = rep(0, nrow(L)), L=diag(2)){
#   return(.Call('_STLM_dmvnorm_chol', PACKAGE='STLM', mu, L))
# }
# 
# rmvt <- function(m = rep(0, nrow(S)), S=diag(2), df=1){
#   return(.Call('_STLM_rmvt_c', PACKAGE='STLM', m, S, df))
# }
# 
# rmvnorm <- function(m = rep(0, nrow(S)), S=diag(2)){
#   return(.Call('_STLM_rmvnorm_c', PACKAGE='STLM', m, S))
# }
# 
# rmvt.chol <- function(m = rep(0, nrow(L)), L=diag(2), df=1){
#   return(.Call('_STLM_rmvt_chol', PACKAGE='STLM', m, L, df))
# }
# 
# rmvnorm.chol <- function(m = rep(0, nrow(L)), L=diag(2)){
#   return(.Call('_STLM_rmvnorm_chol', PACKAGE='STLM', m, L))
# }
# 
# rmvt.chol.R <- function(mu = rep(0, nrow(L)), L=diag(2), df){
#   return(.Call('_STLM_rmvt_chol_R', PACKAGE='STLM', mu, L, df))
# }

rmvt.chol <- function(mu = rep(0, nrow(U)), 
                      U = diag(2), df){
  d <- nrow(U)
  u <- rchisq(1, df)
  
  return(crossprod(U, rnorm(d, 0, 1))/sqrt(u/df) + mu)
}

rmvnorm.chol <- function(mu = rep(0, nrow(U)), U = diag(2)){
  d <- nrow(U)

  return(crossprod(U, rnorm(d, 0, 1)) + mu)
}