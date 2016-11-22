#CHAPTER 18 -- Weibull survival function

library(survival)

#data
t <- c(12.90,2.42,30.35,24.43,13.48,19.49,15.84,8.56,27.15,2.01 ,35.80,49.89,8.66,
17.35,3.24,7.84,9.41,3.91,16.59,37.48 ,24.06,8.08,14.91,6.93,9.46,37.71,20.43,8.06,
23.39,4.68 ,18.11,2.45,14.78,25.88,10.47,9.38,9.62,24.17,20.85,6.26 ,22.09,32.77,
39.17,43.04,6.69,4.92,8.54,10.83,5.56,6.43)

summary(t)

#test -- simulation
#G <- 1.5
#L <- 0.01
#iter <- 10000
#p <- runif(iter)
#t <- ((-log(p))^(1/G))/L

#mean and rate
n <- length(t)
cc <- rep(1,n)
summary(t)
mean(t)
1/mean(t)

f <- survreg(Surv(t,cc)~1,dist="weibull")
summary(f)

s <- f$scale
lambda <- exp(-f$coefficients)
gamma <- 1/f$scale
w.log <- f$loglik[1]
v <- f$var
se <- sqrt(v[2,2])
x2 <- (log(gamma)/se)^2
pvalue <- 1-pchisq(x2,1)
round(cbind(gamma,se,lambda,x2,pvalue),3)

#exponential model
f <- survreg(Surv(t,cc)~1,dist="exponential")
summary(f)
e.log <- f$loglik[1]
lambda.0 <- 1/mean(t)
lambda.1 <- 1/exp(f$coefficients)
cbind(lambda.0,lambda.1)

#likelihood ratio comparison
x2 <- 2*(w.log-e.log)
pvalue <- 1-pchisq(x2,1)
cbind(w.log,e.log,x2,pvalue)

