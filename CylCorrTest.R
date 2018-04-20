library(ggplot2)
library(reshape2)
library(dplyr)

#junk <- matrix(0, nrow=1000, ncol=200)
junk.cone  <- matrix(0, nrow=1000, ncol=100)
junk.cone2 <- matrix(0, nrow=1000, ncol=100)

X <- matrix(seq(1,100, by=1), ncol=1)

# for(i in 1:5000){
#   junk[i, ] <- drawGammaRF(X, 2, 0.5, 1e-3, R=10, sd=i)
# }

for(i in 1:1000){
  junk.cone[i, ]  <- drawGammaRF_1D_Cone_R(X, 2, 0.5, 1e-3, R=10, sd=i+1000)
  junk.cone2[i, ] <- drawGammaRF(X, 2, 0.5, 1e-3, R=10, corrStruct = "Cone", sd=i+2000)
}


# corr.test.1D <- matrix(0, nrow=200, ncol=200)
# for(i in 1:200){
#   for(j in 0:(200-i)){
#     corr.test.1D[i,j] <- cor(junk[,i], junk[, i+j])
#   }
# }
# 
# cyl.test <- melt(corr.test.1D)
# cyl.test <- cyl.test %>% filter(value != 0)
# cyl.test$X <- cyl.test$Var2 * 0.5
# 
# ggplot(data=cyl.test, aes(x=X, y=value)) + 
#   geom_point(alpha=0.1) + 
#   geom_line(data=data.frame(x=X, y=cyl.corr(X, 10, 1)),aes(x=x,y=y), col="red") + 
#   lims(x=c(0, 25)) + 
#   theme_bw()

cone.corr.test.1D   <- matrix(0, nrow=100, ncol=100)
cone.corr.test.1D.2 <- matrix(0, nrow=100, ncol=100)

for(i in 1:100){
  for(j in 0:(100-i)){
    cone.corr.test.1D[i,j]   <- cor(junk.cone[,i], junk.cone[, i+j])
    cone.corr.test.1D.2[i,j] <- cor(junk.cone2[,i], junk.cone2[, i+j])
  }
}

test.cone <- melt(cone.corr.test.1D)
test.cone <- test.cone %>% filter(value != 0)
test.cone$X <- test.cone$Var2 * 0.5

test.cone2 <- melt(cone.corr.test.1D.2)
test.cone2 <- test.cone2 %>% filter(value != 0)
test.cone2$X <- test.cone2$Var2 * 0.5

ggplot(data=test.cone, aes(x=X, y=value)) + 
  geom_point(alpha=0.1) +
  geom_point(data=test.cone2, aes(x=X, y=value), alpha=0.1, col="blue") + 
  geom_line(data=data.frame(x=X, y=cone.corr(X, 5, 1)),aes(x=x,y=y), col="red") + 
  lims(x=c(0, 25)) + 
  theme_bw()



cone.cov.test.1D   <- matrix(0, nrow=200, ncol=200)
cone.cov.test.1D.2 <- matrix(0, nrow=200, ncol=200)

for(i in 1:200){
  for(j in 0:(200-i)){
    cone.cov.test.1D[i,j]   <- cov(junk.cone[,i], junk.cone[, i+j])
    cone.cov.test.1D.2[i,j] <- cov(junk.cone2[,i], junk.cone2[, i+j])
  }
}

test.cone <- melt(cone.cov.test.1D)
test.cone <- test.cone %>% filter(value != 0)
test.cone$X <- test.cone$Var2 * 0.5

test.cone2 <- melt(cone.cov.test.1D.2)
test.cone2 <- test.cone2 %>% filter(value != 0)
test.cone2$X <- test.cone2$Var2 * 0.5

ggplot(data=test.cone, aes(x=X, y=value)) + 
  geom_point(alpha=0.1) +
  geom_point(data=test.cone2, aes(x=X, y=value), alpha=0.1, col="blue") + 
  lims(x=c(0, 25)) + 
  theme_bw()