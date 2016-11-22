#CHAPTER 12 -- Graphics

set.seed(777)
par(mfrow=c(2,1))
par(pty="s")
n <- 1000
cdf <- (1:n)/n
#random exponential distributed data
x <-  rexp(n,1.5)
y <-  rexp(n,1.5)
y0 <- -log(1-cdf)/1.5
#random normal distributed data
#x <-  rnorm(n,2,10)
#y <-  rnorm(n,2,10)
#y0 <- qnorm(cdf,2,10)
plot(sort(x),y0,type="l",xlim=range(x),ylim=range(y))
title("plot quantile functions")
abline(0,1)

#quantile/quantile plot
qqplot(x,y,type="l",xlim=range(x),ylim=range(y))
abline(0,1)
title("qqplot")

#both x and y cumulative probability functions
par(mfrow=c(1,1))
qqplot(x,y,type="l")
lines(sort(x),y0,lty=2)
title("plot cdf-functions")
abline(0,1)



