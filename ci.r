#CHAPTER 2 -- Confidence intervals

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#ci.r
x <- 16
n <- 100
x <- 60
p <- x/n
v <- p*(1-p)/n
ci(p,v)

#function   lower   upper
A <- ci(p,v)[1]
B <- ci(p,v)[2]
c(A,B)
rbind(round(c(100*p,100*c(A,B)),0),
round(c(1-p,1-B,1-A),2),
round(c(1/p,1/B,1/A),2),
round(c(log(p),log(A),log(B)),2),
round(c(log(p/(1-p)),log((A/(1-A))),log((B/(1-B)))),2),
round(c(log(1-p),log(1-B),log(1-A)),2))

#examples of confidence intervals with differing accuracy
#data
p <- c(0.1,0.2,0.5)
n <- 30
x <- n*p

#not adjusted
cbind(p-1.960*sqrt(p*(1-p)/n),p+1.960*sqrt(p*(1-p)/n))
#adjusted
cbind(p-1.960*sqrt(p*(1-p)/n)+1/(2*n),p+1.960*sqrt(p*(1-p)/n)-1/(2*n))
#logistic
v <- ((n+1)*(n+2))/(n*(x+1)*(n-x+1))
l <- log(p/(1-p))
a <- l+1.960*sqrt(v)
b <- l-1.960*sqrt(v)
cbind(1/(1+exp(-b)),1/(1+exp(-a)))
#exact
binom.test(3,n)$conf.int
binom.test(6,n)$conf.int
binom.test(15,n)$conf.int

#median -- confidence interval (n = 15)
#data (ordered):
data  <- c(3.04,3.70,3.80,3.81,3.84,3.89,3.98,4.03,
4.10,4.25,4.26,4.32,4.57,4.84,5.18)
#data <- sort(rnorm(3000,10,5)) #test case
median(data)
n <- length(data)
P <- pbinom(0:n,n,0.5)[-n]
r.data <- rank(data)
round(rbind(P,data,r.data),3)

#bounds and confidence intervals
L <- max(which(pbinom(0:n,n,0.5)<=0.025))
U <- min(which(pbinom(0:n,n,0.5)>=0.975))
cbind(L,U,data[L],data[U])

#95% confidence interval based on the estimated mean
ci(mean(data),var(data)/n)

#large sample (ordered) -- confidence interval (n = 30)
#data (ordered):
x <- c(1.88,2.48,3.02,4.05,6.87,7.82,7.83,8.40,
8.41,8.48,8.78,9.18,9.86,9.90,10.22,10.81,
10.91,11.71,11.96,12.00,12.04,12.49,12.68,13.86,
13.89,14.44,15.09,15.75,16.55,19.25)
#x <- sort(rnorm(1000,10,2)) #test case
median(x)
n <- length(x)
L <- round((n+1)/2-1.960*sqrt(n)/2)
U <- round((n+1)/2+1.960*sqrt(n)/2)
cbind(L,U,x[L],x[U])
ci(mean(x),var(x)/n)



