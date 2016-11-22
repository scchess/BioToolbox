#CHAPTER 22 -- ROC-curve

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#data: test (t) versus standard (s) -- carotid disease
t <- c(20.0,30.0,40.0,45.0,50.0,50.0,55.0,55.0,57.5,
66.7,67.5,67.5 ,73.6,75.0,75.0,75.0,77.5,77.5,80.0,
80.0,80.0,80.0,82.5,82.5,85.0,87.5,87.5,92.5,85.0,85.0)

s <- c(39.5,55.9,59.5,60,63,64.9,66.8,67.2,68.0,70.6,
72.4,73.6 ,75.3,76.3,77.4,77.8,78.2,78.4,78.5,78.6,
80.8,85.6,86.3,87.4,88.5,88.6,88.7,91.5,87.4,88.3)

summary(t)
summary(s)

#probabilities
cbind(1-pnorm(20,18,3),1-pnorm(20,22.5,3))

#"optimum" roc-values -- circle
ns <- 30
nt <- 30
tbar <- mean(t)
sbar <- mean(s)
c <- seq(20,150,0.05)
V <- ((ns-1)*var(s)+(nt-1)*var(t))/(ns+nt-2)
sen <- 1-pnorm(c,sbar,sqrt(V))
fpf <- 1-pnorm(c,tbar,sqrt(V))
plot(fpf,sen,type="l",xlab="false positive",ylab="sensitivity")
abline(0,1)
c0 <- (sbar+tbar)/2
points(cbind(1-pnorm(c0,tbar,sqrt(V)),1-pnorm(c0,sbar,sqrt(V))))
title("ROC curve")

#Parametric ROC analysis
ns <- length(s)
nt <- length(t)
tbar <- mean(t)
sbar <- mean(s)
V <- ((ns-1)*var(s)+(nt-1)*var(t))/(ns+nt-2)
R <- (tbar-sbar)/sqrt(2*V)
auc <- 1-pnorm(R)
z <- (tbar-sbar)/sqrt(V*(1/ns+1/nt))
pvalue <- pnorm(z)
round(cbind(tbar,sbar,V,R,auc,z,pvalue),3)

#95% confidence intervals for R and auc
v.r <- .5*(1/nt+1/ns)
ci(R,v.r) #confidence interval from R
1-pnorm(rev(ci(R,v.r))) #confidence interval from auc

#note:  t-test of difference of mean values tbar and sbar
t.test(t,s,var.equal=TRUE)

#Nonparametric ROC analysis
#artificial data
x <- c(11.49,10.61,11.51,10.60,12.64,11.62,11.20,12.11,10.79,10.62,10.70,11.05,8.12,9.97,12.31,10.97,10.96,9.46,10.67,10.50)
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
n <- length(y)
nt <- sum(y)
ns <- n-nt
U <- wilcox.test(x[y==0],x[y==1])$statistic
W <- U+ns*(ns+1)/2
P <- U/(ns*nt)
wbar <- sum(rank(x)[y==0])/ns
auc <- (wbar-.5*(ns+1))/nt
v <- (1/12)*(1/ns+1/nt)
z <- (auc-0.5)/sqrt(v)
pvalue <- 1-pnorm(z)
round(cbind(U,W,wbar,P,auc,z,pvalue),3)

#data from illustration
s <- c(2.8,5.2,5.6,6.5,8.8,9.9)
t <- c(3.4,5,6,6.1,7,7.3,7.7,8.4,9,9.2)
n <- length(t)+length(s)
nt <- length(t)
ns <- n-nt
cbind(n,ns,nt)
S <- sqrt(((nt-1)*var(t)+(ns-1)*var(s))/(nt+ns-2))
R <- (mean(s)-mean(t))/(S*sqrt(2))
auc <- 1-pnorm(R)
round(cbind(mean(t),mean(s),S,R,auc),3)

#roc plot -- nonparametric
#Table 8.0
ns <- 6
nt <- 10
n <- nt+ns
id <- c(1,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1)
X <- Y <- NULL
x <- y <- 0
for(i in 1:n) {
x <- ifelse(id[i]==1,x+1/ns,x)
X[i] <- 1-x
y <- ifelse(id[i]==0,y+1/nt,y)
Y[i] <- 1-y
}
round(cbind(id,X,Y),2)
plot(c(1,X),c(1,Y),type="s",xlim=c(0,1),ylim=c(0,1),
xlab="false positive fraction",ylab="true positive fraction")
points(c(1,X),c(1,Y),pch=20)
abline(0,1)
title("ROC PLOT")

#a few choices for plot symbols
plot(1:20,rep(5,20),pch=1:20,axes=FALSE,xlab="",ylab="")
text(1:20,rep(4.8,20),1:20)
text(10,5.3,"a few choices for plot symbols")


