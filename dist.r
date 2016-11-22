
#CHAPTER 1 -- Distributions

#probabilities
#normal distribution
z <- seq(-2,2,.5)
round(rbind(z,pnorm(z),1-pnorm(z)),3)

#t-distribution
p <- c(0.990,0.975,0.950,0.900,0.800)
df <- c(2,10,20,40,60)
df
qvalue <- function(p,df) {qt(p,df)}
round(cbind(p,outer(p,df,qvalue)),3)

#chi-square distribution
df <- c(1,2,10,30,50,100)
df
qvalue <- function(p,df) {qchisq(p,df)}
round(cbind(p,outer(p,df,qvalue)),3)

#f-distribution
P <- matrix(0,5,16)
df1 <- c(rep(c(1,2,10,30),c(4,4,4,4)))
df2 <- c(rep(c(10,20,30,60),4))
qvalue <- function(p,df1,df2) {qf(p,df1,df2)}
for(i in 1:5) {
P[i,] <- qvalue(p[i],df1,df2)
}
round(rbind(c(0,df1),c(0,df2),cbind(p,P)),3)

#relationship 1
a <- c(.1,.05,.025,.01)
z2 <- qnorm(1-a/2)^2
Z2 <- qchisq(1-a,1)
cbind(a,z2,Z2)

#relationship 2
df <- c(10,20,50,100)
F <- qf(1-a,1,df)
T  <- qt(1-a/2,df)
cbind(a,df,F,T^2)

#relationship 3
df0 <- c(5,10,60,100)
F <- qf(a,df0,df)
F0 <- qf(1-a,df,df0)
cbind(a,df0,df,F,1/F0)

#relationship 4
df <- c(10,60,150,200)
x2 <- df+qnorm(1-a)*sqrt(2*df)
X2 <- qchisq(1-a,df)
cbind(a,df,x2,X2)

#relationship 5
F <- qf(1-a,df0,df)*df0
X2 <- qchisq(1-a,df0)
cbind(a,df0,df,F,X2)
