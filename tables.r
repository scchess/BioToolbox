#CHAPTER 10 -- Analysis of tables

#birthday data
x <- 148
n <- 348
((x/n-0.5)/sqrt(.25/n))^2
prop.test(x,n,correct=FALSE)
chisq.test(c(x,n-x))$statistic

#birthday -- 12 months
x <- c(25,23,25,25,27,23,42,40,32,35,30,21)
p <- x/n
P <- rep(1/12,12)
round(rbind(x,p,P),3)
G2 <- 2*sum(x*log(p/P))
cbind(G2,1-pchisq(G2,11))
chisq.test(x,p=rep(1/12,12))

#loglikelihood
x <- c(16,12,20,28,24)
p0 <- x/sum(x)
L0 <- sum(x*log(p0))
p1 <- c(0.2,0.2,0.2,0.2,0.2)
L1 <- sum(x*log(p1))
p2 <- c(0.14,0.17,0.20,0.23,0.26)
L2 <- sum(x*log(p2))
p3 <- c(0.14,0.14,0.20,0.26,0.26)
L3 <- sum(x*log(p3))

#loglikelihood test statistics
cbind(L0,L1,L2,L3)
#differences
cbind(2*(L0-L1),2*(L0-L2),2*(L0-L3))

#vitamin C
x <- c(42,27,87,48,29,4)
n <- sum(x)
m <- matrix(x,2,3)
m
e <- outer(apply(m,1,sum),apply(m,2,sum),"*")/n
e
p <- as.vector(m/n)
P <- as.vector(e/n)
round(rbind(p,P),3)
sum(x*log(p))
sum(x*log(P))
G2 <-2*sum(x*log(p/P))
cbind(G2,2*(1-pchisq(G2,2)))
chisq.test(m)


#mean polish -- parity by education
x <- c(26.35,25.36,22.72,25.57,27.49,26.43,26.92,
28.10,23.58,28.17,26.84,26.10)
m <- t(matrix(x,3,4))
m
#model
rr <- factor(rep(1:4,3))
cc <- factor(sort(rep(1:3,4)))
M0 <- as.vector(m)
cbind(M0,rr,cc)
f <- glm(M0~rr+cc)
summary(f)
M <- matrix(fitted(f),4,3)
M
#mean polish
mtemp <- sweep(m,1,apply(m,1,mean),"-")
round(mtemp,2)
mm <- sweep(mtemp,2,apply(mtemp,2,mean),"-")
round(mm,2)
m.add <- m-mm
m.add

#parity by maternal age 
x <- c(3.281,3.305,3.291,3.258,3.225,3.202,3.278,
3.354,3.371,3.356,3.322,3.303,3.280,3.360,
3.397,3.390,3.374,3.335,3.220,3.345,3.405,
3.422,3.412,3.392,3.202,3.332,3.399,3.434,
3.435,3.417,3.333,3.315,3.398,3.443,3.467,3.482)
m <- t(matrix(x,6,6))
m

#nonparametric -- mean polish
mtemp <- sweep(m,1,apply(m,1,mean),"-")
mm <- sweep(mtemp,2,apply(mtemp,2,mean),"-")
#standardized
v <- sum(mm^2)/25
cbind(v,sqrt(v))
round(mm/sqrt(v),3)
ifelse(mm<0,"-","+")

#matched pairs
b <- 225
c <- 132
n <- b+c
X2 <- (b-c)^2/n
X2
prop.test(b,n,correct=FALSE)

#matched -- 5 by 5 table
x <- c(11,15,8,8,6,20,20,12,9,16,12,11,12,8,12,
11,14,13,10,38,28,37,33,46,92)
d <- matrix(x,5,5)
d
#expected
e <- (d+t(d))/2
e
#analysis
X2 <- sum((d-e)^2/e)
X2
#alternative version
x2 <- sum((d-t(d))^2/(d+t(d)))/2
x2
pvalue <- 1-pchisq(X2,10)
cbind(x2,pvalue)

#smoking and chd
#data 
x <- c(29,155,21,76,7,45,12,43)
m <- matrix(x,2,4)
m
#expected values -- independence
e <- outer(apply(m,1,sum),apply(m,2,sum),"*")/sum(m)
e
#model
rr <- factor(c(1,2,1,2,1,2,1,2))
cc <- factor(sort(c(1,2,3,4,1,2,3,4)))
cbind(x,rr,cc)
f <- glm(x~rr+cc,family=poisson)
summary(f)
p <- as.vector(m/sum(m))
P <- as.vector(fitted(f)/sum(m))
n <- as.vector(m)
round(rbind(n,p,P),3)
G2 <- 2*sum(n*log(p/P))
cbind(G2,1-pchisq(G2,3))
#expected values -- model and chisquare
matrix(fitted(f),2,4)

#mean polish
m <- log(m)
mtemp <- sweep(m,1,apply(m,1,median),"-")
mm <- sweep(mtemp,2,apply(mtemp,2,median),"-")
m-mm
exp(m-mm)
#note 
chisq.test(exp(m-mm))$statistic  #independent?

