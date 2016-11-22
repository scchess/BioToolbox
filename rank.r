#CHAPTER 15 -- Nonparametric

#data:
before <- c(131,154,45,89,38,104,127,129,148,181 ,122,
120,66,122,115,113,163,67,156,191)
after  <- c(111,13,16,35,23,115,87,120,122,228 ,118,
121,20,153,139,103,124,37,163,246)

#description
xbar0 <- mean(before)
xbar1 <- mean(after)
n <- length(before)
r <- cor(before,after)
round(cbind(n,xbar0,xbar1,xbar0-xbar1,r),3)
t.test(before-after,correct=F)
mean(before-after)/sqrt(var(before-after)/n)

#Signed rank test
#exact
round(1-pbinom(12,20,0.5),3)
binom.test(13,20,alternative="greater")

#approximate
z <- (0.65-0.5-.5*(1/n))/(0.5*sqrt(1/n))
cbind(z,1-pnorm(z))

#statistical power -- signed rank test
p0 <- 0.5
p <- c(0.5,0.6,0.7,0.8,0.9)
for(n in c(20,30,40,50)){
c <- p0+1.645*0.5*sqrt(1/n)
z <- (c-p)/(0.5*sqrt(1/n))
print(round(c,3))
print(round(1-pnorm(z),2))
}

#Wilcox signed test
wilcox.test(c(17,3,4,5),c(10,8,6,1),paired=T,alternative="greater",exact=T)

#smoking intervention trial data
wilcox.test(before,after,paired=T,alternative="greater",exact=T)
n <- 20
ex <- n*(n+1)/4
v <- n*(n+1)*(2*n+1)/24
W <- wilcox.test(before,after,paired=T,alternative="greater",exact=T)$statistic
z <- (W-ex)/sqrt(v)
x2 <- z^2
pvalue <- (1-pchisq(x2,1))/2
cbind(n,ex,v,W,z,x2,pvalue)

#Kruskal-Wallis
x <- c(2.24,1.59,0.82,9.82,8.55,7.04,5.40,3.73,5.44,
6.58,7.43,8.14,8.76,2.88,4.51,3.05,1.69,0.57,
0.00,7.95,8.13,8.27,9.62,8.55,7.41,6.23,5.04,
3.84,2.76,4.47,11.91,11.52,11.08,10.50,9.75,
2.08,1.09,0.30,0.01,0.58,2.10,3.99,5.62,6.79)

n <- c(14,8,8,14)
id <- rep(1:4,n)
tapply(x,id,length)
tapply(x,id,mean)
tapply(x,id,var)
cbind(sum(n),mean(x),var(x))
summary(aov(x~factor(id)))

#nonparametric
kruskal.test(x,id)
ranks <- rank(x)
tapply(ranks,id,length)
rbar <- tapply(ranks,id,mean)
rbar
tapply(ranks,id,var)
summary(aov(ranks~factor(id)))

N <- length(ranks)
B <- sum(n*(rbar-(N+1)/2)^2)
T <- N*(N+1)*(N-1)/12
X2 <- (N-1)*B/T
pvalue <- 1-pchisq(X2,3)
cbind(T,B,X2,pvalue)

#Three-group regression

xx <- c(133,196,175,205,165,145,198,120,151,180,155,
210,191,146,164,170,157,162,160,135,148)
yy <- c(108,120,128,130,120,104,136,134,110,137,120,
154,125,146,130,148,126,136,127,120,128)

#nonparametric
x <- sort(xx)
m.x <- c(x[4],x[11],x[18])
y <- sort(yy)
m.y <- c(y[4],y[11],y[18])
rbind(m.x,m.y)
B <- (m.y[3]-m.y[1])/(m.x[3]-m.x[1])
A <- sum(m.y)/3-B*sum(m.x)/3
#estimated slope and intercept
cbind(A,B)
#parametric
summary(lm(yy~xx))

#Tukey quick test
#example 1
A <- c(0.2,1.5,2.1,3.0,4.0,4.2,4.4,4.8,5.3,6.7)
B <- c(3.1,3.4,3.7,4.1,5.3,5.8,6.3,6.6,7.8,7.9)
#A contains the minimum and B contains the maximum
t1 <- sum(ifelse(A<min(B),1,0))
t2 <- sum(ifelse(B>max(A),1,0))
t <- t1+t2
cbind(t1,t2,t)

#example 2
B <- c(132,140,149,155,172,179,181,182,185,189,220,
220,225,247,254,255)
A <- c(176,177,177,180,192,202,209,210,211,223,227,
237,238,239,242,271,274,325)

ifelse(B<min(A),1,0)
t1 <- sum(ifelse(B<min(A),1,0))
ifelse(A<max(B),0,1)
t2 <- sum(ifelse(A>max(B),1,0))
t <- t1+t2
cbind(t1,t2,t)

#Friedman test
x <- c(2.3,4.3,5.3,1.6,1.9,6.9,3.5,1.2,0.2)
m <- matrix(x,3,3)
m
friedman.test(m)
M <- t(apply(m,1,rank))
M
rbar <- apply(M,2,mean)
s <- apply(M,2,sum)
rbind(s,rbar)
S <- 3*sum((rbar-mean(rbar))^2)
S

#age by parity -- 6 by 6 table
x <- c(3.281,3.278,3.280,3.220,3.202,3.333,3.305,
3.354,3.360,3.345,3.332,3.315,3.291,3.371,3.397,
3.405,3.399,3.398,3.258,3.356,3.390,3.422,3.434,
3.443,3.225,3.322,3.374,3.412,3.435,3.467,3.202,
3.303,3.335,3.392,3.417,3.482)

m <- t(matrix(x,6,6))
m
M <- t(apply(m,1,rank))
M
friedman.test(m)

rbar <- apply(M,2,mean)
sums <- apply(M,2,sum)
rbind(sums,rbar)
M <- apply(m,1,rank)
R <- apply(M,1,mean)
R
r <- c <- 6
B <- r*sum((R-mean(R))^2)
T <- r*(c-1)*c*(c+1)/12
S <- r*(c-1)*B/T
pvalue <- 1-pchisq(S,5)
cbind(B,T,S,pvalue)

#k by 2 table
#binomial test
r <- 10
c <- 2
prop.test(c,r,correct=FALSE)

#binomial
c1 <- c(1,1,1,2,1,1,2,1,1,1)
c2 <- c(2,2,2,1,2,2,1,2,2,2)
p <- sum(ifelse(c1==2,1,0))/r
x2 <- ((p-.5)/sqrt(1/(4*r)))^2
x2
friedman.test(cbind(c1,c2))

