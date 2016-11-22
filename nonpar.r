#CHAPTER 17 -- Nonparimetric survival analysis

library(survival)

#data
#survival probability distribution (no censoring)
t <- c(3,9,12,15,18,35,48,51,55,68)
n <- length(t)
cc <- rep(1,n)
f <- survfit(Surv(t,cc)~1)
summary(f)
plot(f,xlab="survival time",ylab="survival probability")

#survival probability distribution (censoring)
#data
t <- c(3,9,12,15,18,35,48,51,55,68)
n <- length(t)
cc <- c(1,0,1,1,0,1,0,1,1,1)
f <- survfit(Surv(t,cc)~1)
summary(f)
plot(f,xlab="survival time",ylab="survival probability")

#variance
#complete data
n <- c(10,9,8,7,6,5,4,3,2)
t <- c(6,3,3,3,17,13,3,4,13)

#censored data
n <- c(10,8,7,5,3,2)
t <- c(9,3,20,16,4,13)

#estimation
P <- cumprod(1-1/n)
a <- P*t
nn <- n*(n-1)
A <- rev(cumsum(rev(a)))
AA <- A^2/nn
d <- length(n)+1
v <- (d/(d-1))*sum(AA)
cbind(P,t,a,nn,A,AA)
tbar <- sum(c(1,P)*c(3,t))
cbind(tbar,v)

#complete data: calculated directly from the data
s <- c(3,9,12,15,18,35,48,51,55,68)
mean(s)
var(s)/10

#Cumulative hazard function
t <- c(3,12,15,35,51,55,68)
n <- c(10,8,7,5,3,2,1)
p <- 1-1/n
H <- cumsum(-log(p))
round(cbind(t,n,p,H),3)
plot(t,H,type="s",ylim=c(0,2),xlim=c(0,68),xlab="Time",ylab="Cumulative hazard function")

#log-rank analysis
#data -- SFMHS -- smokers (n = 8) and nonsmokers (n = 17)
smk <- c(2,8,13,24,25,38,43,48)
c.smk <- c(1,1,1,0,1,1,1,1)  #censored = 0
smk0 <- c(12,18,21,34,40,46,30,33,39,42,44,50,56,58,61,66,77)
c.smk0 <- c(1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,0,1)
t <- c(smk,smk0)
cc <- c(c.smk,c.smk0)
status <- rep(0:1,c(length(smk),length(smk0)))

#log-rank analysis
f <- survfit(Surv(t,cc)~status)
summary(f)
f0 <- survdiff(Surv(t,cc)~status)
f0
cbind(f0$exp[1],f0$var[1,1],f0$chisq)


