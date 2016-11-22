#CHAPTER 11 -- Bootstrap estimation

#smoking intervention trial
#data
set.seed(777)
before <- c(131,154,45,89,38,104,127,129,148,181,122,
120,66,122,115,113,163,67,156,191)
after <- c(111,13,16,35,23,115,87,120,122,228,118,121,
20,153,139,103,124,37,163,246)
P <- 100*(before-after)/before
summary(P)
round(rbind(before,after,P),1)

t.test(before,after,paired=TRUE,var.equal=TRUE)

#bootstrap estimation -- intervention trial 
iter <- 2000
PBAR <- NULL
P <- 100*(before-after)/before
for(i in 1:iter) {
ptemp <- sample(P,length(P),replace=TRUE)
PBAR[i] <- mean(ptemp)
}
pbar <- mean(PBAR)
vbar <- var(PBAR)
sd <- sqrt(vbar)
cbind(pbar,vbar,sd)
hist(PBAR,50)

#bootstrap test statistics
pvalue.0 <- sum(ifelse(PBAR<0,1,0))/iter
pvalue.1 <- pnorm((0-pbar)/sd)
cbind(pvalue.0,pvalue.1)

#bootstrap confidence intervals
cbind(pbar-1.96*sd,pbar+1.96*sd)
cbind(sort(PBAR)[.025*iter],sort(PBAR)[.975*iter])

#example bootstrap for tabled data
iter <- 2000
p <- NULL
d <- rep(1:3,c(10,5,5))
n <- length(d)
for(i in 1:iter) {
D <- sample(d,n,replace=TRUE)
P <- table(factor(D,level=1:3))
p[i] <- (P[1]-2*P[2]+P[3])/n
}
cbind(mean(p),var(p))

#Kappa statistic --- a bootstrap estimate
a <- 105
b <- 24
c <- 18
d <- 165
n <- a+b+c+d
data <- rep(1:4,c(a,b,c,d))
table(data)/n
iter <- 2000
K <- NULL
for (i in 1:iter ) {
temp <- sample(data,n,replace=TRUE)
tab <- table(temp)
P <- (tab[1]+tab[4])/n 
q1 <- (tab[1]+tab[2])/n
q2 <- 1-q1
p1 <- (tab[1]+tab[3])/n
p2 <- 1-p1
p <- p1*q1+p2*q2
K[i] <- (P-p)/(1-p)
}
m <- mean(K)
v <- var(K)
cbind(m,v)
#parametric
lower <- m-1.96*sqrt(v)
upper <- m+1.96*sqrt(v)
cbind(lower,m,upper)
#nonparametric
LOWER <- sort(K)[iter*0.025]
UPPER <- sort(K)[iter*0.975]
cbind(LOWER,m,UPPER)
hist(K,50)


