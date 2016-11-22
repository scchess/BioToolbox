#CHAPTER 6 Two by two tables
#twobytwo.6

#Wilcoxon test
temp <- c(344,0,246,0,224,1,242,0,252,1,233,1,224,0,239,1,252,0,202,1,291,1,212,0,239,1,153,0,218,1,312,1,188,0,254,1,183,0,202,0,185,0,250,0,169,0,234,1,212,1,250,1,197,1,226,0,137,0,325,1,263,0,148,0,175,0,181,1,194,0,246,1,268,1,276,1,248,1,213,0)

d <- temp[seq(1,length(temp),2)]
id <- temp[seq(2,length(temp),2)]
ca <- d[id==1]
cb <- d[id==0]
na <- length(ca)
nb <- length(cb)
cbind(mean(ca),mean(cb),var(ca),var(cb),mean(ca)-mean(cb))

U <- wilcox.test(ca,cb)$statistic
wilcox.test(cb,ca,correct=FALSE)

W <- U+na*(na+1)/2
P <- U/(na*nb)
cbind(U,W,P)

#chd and smoking
r1 <- c(25,24,20,46,45,50,47)
r0 <- c(417,426,327,493,350,398,485)
x <- 0:6

#cancer and xrays
r1 <- c(7332,287,199,96,59,65)
r0 <- c(7673,239,154,65,28,29)
x <- 0:5

#analysis
n1 <- sum(r1)
n0 <- sum(r0)
n <- n1+n0
r <- r1+r0
xbar1 <- sum(x*r1)/sum(r1)
xbar0 <- sum(x*r0)/sum(r0)
xbar <- (n1*xbar1+n0*xbar0)/n
cbind(xbar1,xbar0,xbar,xbar1-xbar0)
sxx <- sum(r*(x-xbar)^2)
syy <- n1*n0/n
sxy <- (xbar1-xbar0)*syy
b <- sxy/sxx
P <- n1/n
a <- P-b*xbar
cbind(a,b,P)
p <- r1/r
pp <- a+b*x
rbind(p,pp)

X2 <- sum(r*(p-P)^2)/(P*(1-P))
XL <- sum(r*(pp-P)^2)/(P*(1-P))
XNL <- sum(r*(p-pp)^2)/(P*(1-P))
cbind(XNL,XL,X2,XNL+XL)

#analysis of rows ---difference in mean values
s <- sqrt((sxx/(n0+n1))*(1/n0+1/n1))
z <- (xbar1-xbar0)/s
cbind(s,z,z^2)

#small example table
x <- c(9,32,94,119,53,74,60,82)
m <- matrix(x,2,4)
m
chisq.test(m)
