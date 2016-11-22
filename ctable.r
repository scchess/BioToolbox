#CHAPTER 8 -- Contingency tables

#data
wilcox.test(c(4,35,21,28,66),c(10,42,71,77,90))
U <- wilcox.test(c(4,35,21,28,66),c(10,42,71,77,90))$statistic
n1 <- 5
W <- U+n1*(n1+1)/2
P <- U/(n1*n1)
cbind(n1,U,W,P)

#A/B behavior and cholesterol
x <- c(344,0,246,0,224,1,242,0,252,1,233,1,224,0,
239,1,252,0,202,1,291,1,212,0,239,1,153,0,218,1,
312,1,188,0,254,1,183,0,202,0,185,0,250,0,169,0,
234,1,212,1,250,1,197,1,226,0,137,0,325,1,263,0,
148,0,175,0,181,1,194,0,246,1,268,1,276,1,248,1,213,0)

matrix(x,2)

chol <-x[seq(1,length(x),2)]
type <-x[seq(2,length(x),2)]
t.test(chol[type==0],chol[type==1])
t.test(chol[type==0],chol[type==1],var.equal=TRUE)
wilcox.test(chol[type==1],chol[type==0])
wilcox.test(chol[type==1],chol[type==0],correct=FALSE)
U <- wilcox.test(chol[type==1],chol[type==0])$statistic
n <- length(type)
n1 <- sum(type)
n2 <- n-n1
W <- U+n1*(n1+1)/2
P <- U/(n1*n2)
v <- (n+1)/(12*n1*n2)
z <- (P-0.5)/sqrt(v)
pvalue <- 2*(1-pnorm(z))
cbind(W,U,P,v,z,pvalue)
#check
sum(rank(chol)[type==0])
sum(rank(chol)[type==1])
U+n1*(n1+1)/2
sum(rank(chol)[type==0])+sum(rank(chol)[type==1])
n*(n+1)/2


#analysis 2 by k table -- Mann/Whitney
x <- c(7,25,7,34,8,32,20,32)
m <- matrix(x,2,4)
m
y <- rep(1:4,m[1,])
n <- rep(1:4,m[2,])
M <- outer(y,n,"-")
U <- sum(table(M[M>0]))+length(M[M==0])/2
N <- length(M)
P <- U/N
cbind(U,N,P)

#example data
n1 <- c(18,16,24,24,52)
n2 <- c(20,24,60,35,50)
x <- 0:4

#body weight data (data set 1)
#n1 <- c(25,24,20,46,45,50,47)
#n2 <- c(417,426,327,493,350,398,485)
#x <- 0:6

#x-ray data (data set 2)
n1 <- c(7332,287,199,96,59,65)
n2 <- c(7673,239,154,65,28,29)
x <- 0:5


#analysis
n <- n1+n2
xbar1 <- sum(x*n1)/sum(n1)
xbar2 <- sum(x*n2)/sum(n2)
xbar <- (sum(n1)*xbar1+sum(n2)*xbar2)/sum(n)
cbind(xbar1,xbar2,xbar)
sxx <- sum(n*(x-xbar)^2)
N1 <- sum(n1)
N2 <- sum(n2)
v <- (sxx/(N1+N2))*(1/N1+1/N2)
z <- (xbar1-xbar2)/sqrt(v)
cbind(z^2,1-pchisq(z^2,1))
syy <- sum(n1)*sum(n2)/sum(n)
sxy <- (xbar1-xbar2)*syy
cbind(syy,sxx,sxy)
b <- sxy/sxx
P <- sum(n1)/sum(n)
a <- P-b*xbar
cbind(a,b,P)
p <- n1/n
pp <- a+b*x
round(rbind(p,pp),3)
X2 <- sum(n*(p-P)^2)/(P*(1-P))
XL <- sum(n*(pp-P)^2)/(P*(1-P))
XNL <- sum(n*(p-pp)^2)/(P*(1-P))
cbind(X2,XL,XNL)


#Comparison of mean values  
#example data
d1 <- c(18,16,24,24,52)
d2 <- c(20,24,60,35,50)
x <- 0:4

#body weight data (data set 1)
d1 <- c(25,24,20,46,45,50,47)
d2 <- c(417,426,327,493,350,398,485)
x <- 0:6

#maternal x-ray data (data set 2)
d1 <- c(7332,287,199,96,59,65)
d2 <- c(7673,239,154,65,28,29)
x <- 0:5

#analysis of mean values
n <- d1+d2
xbar1 <- sum(x*d1)/sum(d1)
xbar2 <- sum(x*d2)/sum(d2)
xbar <- sum(x*(d1+d2))/sum(n)
cbind(xbar1,xbar2,xbar)
sxx <- sum(n*(x-xbar)^2)
v <- (sxx/sum(n))*(1/sum(d1)+1/sum(d2))
z <- (xbar1-xbar2)/sqrt(v)
cbind(z,z^2)

