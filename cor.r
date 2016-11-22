#CHAPTER 5 -- Correlation

#data
x <- c(144,154,142,152,250,204,164,259,190,175,167,186,153,166,180,251,168,225,196,160)
y <- c(77,79,78,74,86,89,84,87,76,90,91,92,73,85,80,88,83,81,71,82)

cbind(mean(x),mean(y),sd(x),sd(y),cov(x,y))
cor.test(x,y)
cor.test(x,y)$estimate
cor.test(x,y,method="spearman")$estimate
cor.test(x,y,method="kendall")$estimate

#estimated x/y-slope
f <- lm(y~x)
summary(f)
summary(f)$coefficient[2,]
cor(x,y)*sd(y)/sd(x)

#ranked data
X <- rank(x)
Y <- rank(y)
cbind(mean(X),mean(Y),sd(X),sd(Y),cov(X,Y))
cov(X,Y)/(sd(X)*sd(Y))
cor.test(x,y,method="spearman")

#point serial correlation coefficient
x <- c(225,177,181,132,255,182,155,140,149,325,223,271,
238,189,140,247,220,176,185,202,227,192,179,
237,177,239,210,220,274,225)
y <- c(1,1,0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,1,0,1,1,1,0,1,1,1,1,0,1,0)
cbind(mean(x),mean(y),sd(x),sd(y),cov(x,y))
x0 <- x[y==0]
y0 <- x[y==1]
cbind(mean(y0),mean(x0),
sum((x-mean(x))^2),
sum((y-mean(y))^2),
sum((x-mean(x))*(y-mean(y))))
cor(x,y)

#short-cut expressions for ssy and sxy
n0 <- sum(y)
n1 <-length(y)-n0
n <- n0+n1
cbind(n0*n1/n,sum((y-mean(y))^2))
cbind(n1*n0*(mean(x[y==1])-mean(x[y==0]))/n,sum((x-mean(x))*(y-mean(y))))

#equivalent tests
cor.test(x,y)
t.test(x[y==1],x[y==0],var.equal=TRUE)
cor.test(x,y,continuity=TRUE)$statistic
t.test(x[y==1],x[y==0],var.equal=TRUE)$statistic

#point biserial correlation
r.pb <- sqrt(n0*n1/n)*(mean(x[y==1])-mean(x[y==0]))/sqrt(sum((x-mean(x))^2))
cbind(r.pb,cor.test(x,y)$estimate)

#gamma-coefficient

#data 1
#artifical example data (n = 4)
x <- c(1,6,5,3)
y <- c(2,5,3,6)

#data 2
#blood pressure data (n = 20)
x <- c(144,154,142,152,250,204,164,259,190,175,167,186,153,166,180,251,168,225,196,160)
y <- c(77,79,78,74,86,89,84,87,76,90,91,92,73,85,80,88,83,81,71,82)

S <- sign(outer(x,x,"-")*outer(y,y,"-"))
C <- length(S[S>0])
D <- length(S[S<0])
gamma <- (C-D)/(C+D)
cbind(C,D,gamma)
cor.test(x,y,method="kendall")

#proportional reduction in error criterion (3 by 5 table)
#lambda correlation coefficient

x <- c(78,102,103,124,153,81,149,220,32,221,352,48,343,481,35)
m <- matrix(x,3,5)
m
P0 <- 1-max(apply(m,2,sum))/sum(m)
p.row <- 1-apply(m,1,max)/apply(m,1,sum)
P <- sum(apply(m,1,sum)*p.row)/sum(m) 
lambda <- (P0-P)/P0
round(p.row,3)
round(cbind(P0,P,lambda),3)

