#CHAPTER 4 -- Binomial/Poisson distributions

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#binomial
n <- 10
p0 <- 0.5
x <- 6
rbind(0:10,round(dbinom(0:n,n,p0),3))
rbind(0:10,round(1-pbinom(0:n,n,p0),3))

#three estimates of the same probability -- P(greater than x = 5)
1-pbinom(x,n,p0)
1-pnorm((x-n*p0+.5)/sqrt(n*p0*(1-p0)))  #approximate
1-binom.test(x,n,p=p0,alternative = "less")$p.value  #exact

#note:
n <- 9
choose(n,0:n)
cbind(sum(choose(n,0:n)),2^n)

#Viral prevalence
k <- 60
x <- 18
n <- 20
Q <- x/k
p <- 1-Q^(1/n)
V <- Q*(1-Q)/k
v <- (Q^(1/n)/(n*Q))^2*V
cbind(Q,p,v,V)
ci(Q,V)
rev((1-ci(Q,V)^(1/n)))
ci(p,v)

#randomized response
x <- 127
n <- 200
pi <- 0.3
p <- x/n
P <- (p-(1-pi))/(2*pi-1)
v <- (p*(1-p)/n)/(2*pi-1)^2
cbind(p,P,v)
ci(P,v)

#geometric series
k <- 0:8
p <- 0.4
q <- 1-p
round(p*q^k,3)
#expected value (q/p)
k <- 1:80
cbind(sum(k*p*q^k),q/p)


#Poisson -- spatial distribution
p <- c(dpois(0:5,2),1-ppois(5,2))
e <- 100*p
o <- c(16,27,23,17,12,2,3)
round(rbind(p,e,o),3)
X2 <- sum((o-e)^2/e)
pvalue <- 1-pchisq(X2,6)
cbind(X2,pvalue)

#goodness-of-fit
k <- 0:6 
x <- log(e)+log(factorial(k))
y <- -2+log(2)*k
plot(x,y,type="b",xlab="Log(expected value)",ylab="Log(observed)")
title("observed/theoretical -- Poisson distributions ?")
x <- log(o)+log(factorial(k))
lines(x,y,type="b",pch=16,lty=2)

#truncated Poisson distribution
#data
d <- c(2,2,2,2,1,1,2,2,2,3,2,4,1,2,3,2,5,3,1,1,2,2,1,2,2,2,1,2,2,4,3,2,2,6,1,5,2,
4,1,4,2,6,2,2,2,1,5,4,3,4,1,4,3,3,2,3,1,3,1,2,3,4,3,2,3,3,1,2,5,1,3,1,2,3,3,3,3,4,
3,5,1,3,3,3,3,1,4,1,6,1,2,1,3,2,3,2,4,2,4,2,2,2,1,3,3,3,6,4,3,6,4,2,3,3,3,4)

tab <- table(d)
tab
mean(d)
N <- sum(d)
S <- tab[1]
n <- length(d)
cbind(N,S,n)
lambda <- (N-S)/n
cbind(N,S,n,lambda)

#maximum likelihood estimation
f <- function(x) {
mean(d)-x/(1-exp(-x))
}
uniroot(f,c(1,4))$root

#alternatively
x <- seq(1,4,0.0001)
F <- abs(f(x))
k <- which.min(F)
x[k]
