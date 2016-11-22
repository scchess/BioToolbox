#CHAPTER 14 -- Log-normal distribution

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#example (Figure 1.0) --  P(Y>4.48) = P(X>1.5)
1-plnorm(4.48169,1,sqrt(0.7))
1-plnorm(exp(1.5),1,sqrt(0.7))
1-pnorm(1.5,1,sqrt(0.7))

y <- rlnorm(10^3,0.5,1.0)
hist(y,50,freq=FALSE)
lines(density(y))
hist(log(y),50,freq=FALSE)
lines(density(log(y)))

#simulated lognorm data
y <- c(2.69,1.11,2.75,1.11,8.49,3.07,2.02,5.00,1.34,
1.13,1.22,1.74,0.25,1.59,16.63,4.36,4.33,0.96,3.23,
2.72,0.22,2.07,0.75,5.89,6.94,2.53,0.29,1.61,0.37,0.96)

#more simulated lognormal data -- mean = 0.5 and variance = 1
#n <- 10000
#y <- rlnorm(n,0.5,1)

x <- log(y)
n <- length(x)
m <- mean(x)
v <- var(x)
ci(m,v/n)
exp(ci(m,v/n))

#percentiles
c <- seq(0.2,0.9,0.1)
round(cbind(c,exp(m+qnorm(c)*sqrt(v))),3)

#mean variance | mean, variance, median and mode
round(cbind(m,v,exp(m+0.5*v),
exp(2*m+v)*(exp(v)-1),
exp(m),
exp(m-sqrt(v))),3)
