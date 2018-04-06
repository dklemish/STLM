plot.circle <- function(c, r, Color="blue"){
  x.axis <- seq(c[1]-r, c[1]+r, by=r/100)
  lines(x.axis, c[2] + sqrt(r^2-(x.axis-c[1])^2), col=Color)
  lines(x.axis, c[2] - sqrt(r^2-(x.axis-c[1])^2), col=Color)
}

X <- matrix(c(0.25, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.75), nrow=4)
area <- matrix(runif(10e6, -0.25, 1.25), ncol=2)

valid.new.pt <- rep(TRUE, 5e6)
valid.new.pt <- valid.new.pt & (((area[,1]-X[1,1])^2 + (area[,2]-X[1,2])^2) < 0.25)
num1 <- sum(valid.new.pt)
sample1 <- area[valid.new.pt,][sample.int(sum(valid.new.pt), 1000),]

valid.new.pt <- rep(TRUE, 5e6)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[1,1])^2 + (area[,2]-X[1,2])^2) < 0.25)
valid.new.pt <- valid.new.pt &  (((area[,1]-X[2,1])^2 + (area[,2]-X[2,2])^2) < 0.25)
num2 <- sum(valid.new.pt)
sample2 <- area[valid.new.pt,][sample.int(sum(valid.new.pt), 1000),]

valid.new.pt <- rep(TRUE, 5e6)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[1,1])^2 + (area[,2]-X[1,2])^2) < 0.25)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[2,1])^2 + (area[,2]-X[2,2])^2) < 0.25)
valid.new.pt <- valid.new.pt &  (((area[,1]-X[3,1])^2 + (area[,2]-X[3,2])^2) < 0.25)
num3 <- sum(valid.new.pt)
sample3 <- area[valid.new.pt,][sample.int(sum(valid.new.pt), 1000),]

valid.new.pt <- rep(TRUE, 5e6)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[1,1])^2 + (area[,2]-X[1,2])^2) < 0.25)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[2,1])^2 + (area[,2]-X[2,2])^2) < 0.25)
valid.new.pt <- valid.new.pt & !(((area[,1]-X[3,1])^2 + (area[,2]-X[3,2])^2) < 0.25)
valid.new.pt <- valid.new.pt &  (((area[,1]-X[4,1])^2 + (area[,2]-X[4,2])^2) < 0.25)
num4 <- sum(valid.new.pt)
sample4 <- area[valid.new.pt,][sample.int(sum(valid.new.pt), 1000),]

plot(sample1, xlim=c(-0.25, 1.25), ylim=c(-0.25, 1.25), pch=16)
lines(sample2, pch=16, col="blue", type="p")
lines(sample3, pch=16, col="red", type="p")
lines(sample4, pch=16, col="green", type="p")

(num1 / 5e6) * 1.5^2
(num2 / 5e6) * 1.5^2
(num3 / 5e6) * 1.5^2
(num4 / 5e6) * 1.5^2
