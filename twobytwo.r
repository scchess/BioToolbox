
#CHAPTER 6 -- Two by two tables

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#typical summary statistics: odds ratio, relative risk and attributable risk
a <- 61
b <- 102
c <- 52
d <- 196
n <- a+b+c+d
n1 <- a+b
n2 <- c+d
m1 <- a+c
m2 <- b+d

or <- (a*d)/(b*c)
rr <- (a/(a+b))/(c/(c+d))
ar <- ((a+c)/n-(c/(c+d)))/((a+c)/n)
cbind(or,rr,ar)

#association measured by probabilities
#rows
prop.test(c(61,52),c(163,248),correct=FALSE)
#columns
prop.test(c(61,102),c(113,298),correct=FALSE)
chisq.test(matrix(c(61,52,102,196),2,2),correct=FALSE)

#Table 7.0 (toolbox text)
p1 <- a/(a+b)
p2 <- c/(c+d)
P1 <- a/(a+c)
P2 <- b/(b+d)
vp <- p1*(1-p1)/n1+p2*(1-p2)/n2
vP <- P1*(1-P1)/m1+P2*(1-P2)/m2
vrr <- rr^2*(1/a-1/n1+1/c-1/n2)
lvrr <- 1/a-1/n1+1/c-1/n2
vor <- or^2*(1/a+1/b+1/c+1/d)
lvor <- 1/a+1/b+1/c+1/d
var0 <- (1-ar)^2*(b+ar*(a+d))/(n*c)
lvar0 <- (b+ar*(a+d))/(n*c)
round(cbind(sqrt(vp),sqrt(vP),sqrt(vrr),sqrt(lvrr),
sqrt(vor),sqrt(lvor),sqrt(var0),sqrt(lvar0)),3)

#Table 8.0 (toolbox text)
p1 <- a/(a+b)
p2 <- c/(c+d)
P1 <- a/(a+c)
P2 <- b/(b+d)
rr <- (a/(a+b))/(c/(c+d))
lrr <- log(rr)
or <- (a/b)/(c/d)
lor <- log(or)
ar <- (a*d-b*c)/((a+c)*(c+d))
lar <- log(1-ar)
round(cbind(p1-p2,P1-P2,rr,lrr,or,lor,ar,lar),3)

#confidence intervals
rbind(round(c(p1-p2,ci(p1-p2,vp)),3),
round(c(P1-P2,ci(P1-P2,vP)),3),
round(c(exp(lrr),exp(ci(lrr,lvrr))),3),
round(c(lrr,ci(lrr,lvrr)),3),
round(c(exp(lor),exp(ci(lor,lvor))),3),
round(c(lor,ci(lor,lvor)),3),
round(c(1-exp(lar),1-exp(rev(ci(lar,lvar0)))),3),
round(c(lar,ci(lar,lvar0)),3))

#corrected
p <- (a+c)/(a+b+c+d)
v <- p*(1-p)*(1/n1+1/n2)
z <- (abs(p2-p1)-.5*(1/n1+1/n2))/sqrt(v)
cbind(p,v,z^2)
prop.test(c(a,c),c(n1,n2),correct=TRUE)$statistic
n*(abs(a*d-b*c)-n/2)^2/((a+b)*(c+d)*(a+c)*(b+d))

#corrected estimate
p <- c(1,4,6,4,1)/16
ex <- sum((0:4)*p)
c0 <- 1:5
p0 <- cumsum(p)  #exact
p1 <- pnorm(c0-0.5-ex)  #approximate
round(cbind(c0-1,p0,p1),3)
chisq.test(matrix(c(a,c,b,d),2,2),correct=FALSE)$statistic
prop.test(c(a,c),c(n1,n2),correct=FALSE)$statistic
n*(abs(a*d-b*c))^2/((a+b)*(c+d)*(a+c)*(b+d))

#Vietnam -- breast cancer
a <- 170
b <- 3222
c <- 126
d <- 2912
n <- a+b+c+d
#relative risk
(a/(a+b))/(c/(c+d))
#odds ratio
(a/b)/(c/d)
#variance of log-odds-ratio
1/a+1/b+1/c+1/d

#Adjustment  #small sample size
a <- 2
b <- 23
c <- 6
d <- 22
n <- a+b+c+d
m <- matrix(c(a,c,b,d),2,2)
chisq.test(m,correct=FALSE)
chisq.test(m)

#Hypergeometric probability distributions

#tea tasting experiment
round(dhyper(0:4,4,4,4),3)
#keno
round(dhyper(0:8,20,60,8),6)
#Fisher's exact test
round(dhyper(0:8,8,20,8),4)
round(1-phyper(0:8,8,20,8),4)

#Adjusted/unadjusted -- Fisher's exact test
a <- 5
b <- 3
c <- 3
d <- 17
n <- a+b+c+d
m <- matrix(c(a,c,b,d),2,2)
chisq.test(m,correct=FALSE)
chisq.test(m)
fisher.test(m)

#twins
#no birth defects (data set 1) 
a <- 18687
b <- 16093
c <- 19188
#birth defects  (data set 2)
a <- 168
b <- 53
c <- 110
n <- a+b+c
p <- (2*a+b)/(2*n)
q <- 1-p
r <- 1-b/(2*p*q*n)
vr <- ((1-r)*(1-2*p*q*(1-r)-(1-4*p*q)*(1-r)^2))/(2*p*q*n)
vp <- (2*p*q*(1+r))/(2*n)
cbind(n,p,r,vr,vp)
ci(r,vr)

#duffy
a <- 8
b <- 72
c <- 409
n <- a+b+c
p <- (2*a+b)/(2*n)
q <- 1-p
chisq.test(c(a,b,c),p=c(p^2,2*p*(1-p),(1-p)^2))
r <- 1-b/(2*p*q*n)
X2 <- n*r^2
pvalue <- 1-pchisq(X2,1)
cbind(p,r,X2,pvalue)  #degrees of freedom = 1)

#incomplete data -- zero in a 2 by 2 table
a <- 33
b <- 44
c <- 14
n <- a+b+c
d <- b*c/a
n+d
(a+b)*(a+c)/a

#false positive
d <- c(0.5,0.1,0.05,0.005,0.0005)
t <- (.95*d+.1*(1-d))
p <- .1*(1-d)/t
round(cbind(d,t,p),3)
