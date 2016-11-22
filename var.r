#CHAPTER 13 -- Variance

#confidence interval
#data
x <- c(10.92,9.20,11.02,9.20,13.48,11.24,10.41,12.22,9.59,9.24) 
mean(x)
var(x)
df <- length(x)-1
lower <- df*var(x)/qchisq(0.975,df)
upper <- df*var(x)/qchisq(0.025,df)
cbind(lower,mean(x),upper)

#test of variance -- Poisson?
#data
x <- c(1,3,1,0,1,0,0,1,0,2,2,0,0,0,1,1,0,1,2,2,0,
1,2,1,1,0,2,2, 1,0,2,1,1,1,1,0,0,1,1,3,2,
0,0,4,3,1,0,1,2,2,1,2,1,1,0,0,1,1,2,2)
tab <- table(x)
tab
X2 <- sum((x-mean(x))^2/mean(x))
pvalue <- 1-pchisq(X2,59)
cbind(mean(x),var(x),X2,pvalue)


#birth weight and birth order --- analysis of variance
wt <- c(2.48,2.92,3.73,4.16,3.80,3.42,3.62,3.82,3.92,3.34,3.26,
3.87,2.92,3.20,4.10,4.06,3.35,3.40,4.48,4.00,3.42,2.76,
2.98,4.58,4.23,4.26,3.66,4.13,2.80,3.62,3.40,3.70,3.10,
3.56,4.08,2.38,3.74 ,3.56,3.62,3.52)
order <- c(0,1,3,1,1,2,0,2,2,1,0,1,1,0,1,0,0,1,0,1,1,0,0,2,3,0,0,3,0,1,1,0,0,0,0,1,0,0,0,1)
#summary tables
tapply(wt,order,length)
tapply(wt,order,mean)
tapply(wt,order,var)
cbind(length(wt),mean(wt),var(wt))
#analysis of variance
summary(aov(wt~factor(order)))

#five groups...homogeneity/heterogeneity
x <- c(1,5,14,24,72)
n <- c(10,25,35,40,90)
p <- x/n
v <- n*p*(1-p)
N <- sum(n)
P <- sum(x)/N
round(cbind(x,n,p,v),2)
cbind(sum(x),sum(n),sum(n*p)/N,N*P*(1-P))
N*P*(1-P)
sum(n*p*(1-p))
sum(n*(p-P)^2)

#tests of variances
#f-test
chol <- c(137,148,153,169,175,181,183,185,188,194,197,
202,202,212,212,213,218,224,224,226,233,234,
239,239,242,246,246,248,250,250,252,252,254,
263,268,276,291,312,325,344)
ab <- c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,0,1,
1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,1,1,1,0)
ranks <- c(1,4,5,8,9,12,13,16,17,20,21,24,25,28,
29,32,33,36,37,40,39,38,35,34,31,30,27,26,23,22,19,18,
15,14,11,10,7,6,3,2)

tapply(chol,ab,length)
tapply(chol,ab,mean)
tapply(chol,ab,var)
var.test(chol[ab==0],chol[ab==1])

#Bartlett test
bartlett.test(chol,ab)
n <- 40
n0 <- 20
n1 <- 20
v.0 <- var(chol[ab==0])
v.1 <- var(chol[ab==1])
v.w <- ((n0-1)*v.0+(n1-1)*v.1)/(n-2)
X2 <- (n-2)*log(v.w)-((n0-1)*log(v.0)+(n1-1)*log(v.1))
#c = correction factor
c <- 1+(1/3)*(2/19-1/38)
X2.c <- X2/c
pvalue <- 1-pchisq(X2.c,1)
round(cbind(v.w,v.0,v.1,X2,X2.c,pvalue),3)

#Levene test
y.a <- abs(chol[ab==1]-mean(chol[ab==1]))
y.b <- abs(chol[ab==0]-mean(chol[ab==0]))
cbind(mean(y.a),mean(y.b))
t.test(y.b,y.a)

#Siegel/Tukey
W <- sum(ranks[ab==1])
n1 <- length(ranks[ab==1])
W0 <- n1*(2*n1+1)/2 
v <- n1*n1*(2*n1+1)/12
cbind(W,W0,v)
X2 <- ((W-W0)/sqrt(v))^2
pvalue <- 1-pchisq(X2,1)
cbind(X2,pvalue)
