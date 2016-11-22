#CHAPTER 3 -- weighted averages

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#least squares versus weighted average
x <- c(1,9,8,6,3,1,3,9)
y <- c(3,8,9,7,2,1,6,4)

#least squares estimation
f <- glm(y~x)
b.lsq <-f$coefficient[2]
#weighted average
b.wt <- sum((x-mean(x))*(y-mean(x)))/sum((x-mean(x))^2)
cbind(b.lsq,b.wt)

#Childhood leukemia  -- confidence interval
data<-c(33,4731047,181,18368181,231,23326046,260,23882120,245,23156490)
deaths <- data[seq(1,length(data),2)]
pop <- data[seq(2,length(data),2)]
rate <- deaths/pop
round(cbind(pop,deaths,rate*10^6),2)
lrate <- log(rate)
lbar <- sum(deaths*lrate)/sum(deaths)
rate <- exp(lbar)*10^6
v <- 1/sum(deaths)
cbind(lbar,ci(lbar,v))
cbind(rate,exp(ci(lbar,v))*10^6)

#weighted.r -- rate ratio
data <- c(33,4731047,238,36243470,181,18368181,1118,142332392,231,
23326046,1253,179720009,260,23882120,1505,186897294,245,23156490,1950,185997993)
matrix(data,4)

#Leukemia data California and US
d <- data[seq(1,length(data),4)]
p <- data[seq(2,length(data),4)]
r <- 10^6*d/p
D <- data[seq(3,length(data),4)]
P <- data[seq(4,length(data),4)]
R <- 10^6*D/P
round(cbind(d,p,r,D,P,R),2)

#confidence interval leukemia California and US
rr <- (d/p)/(D/P)
v <- 1/d+1/D
w <- 1/v
lrr <- sum(w*log(rr))/sum(w)
vrr <- 1/sum(w)
rate <- exp(lrr)
cbind(lrr,vrr,ci(lrr,vrr))
cbind(rate,exp(ci(lrr,vrr)))

#summary odds ratio
a <- c(98,54,11,7)
b <- c(832,227,85,102)
c <- c(169,55,61,90)
d <- c(3520,686,926,1936)
n <- a+b+c+d

#Mantel/Haenszel summary odds ratio
or.mh <- sum(a*d/n)/sum(b*c/n)
or.mh

#weighted average summary odds ratio
or <- a*d/(b*c)
lor <- log(or)
v <- 1/a+1/b+1/c+1/d
w <- 1/v
round(cbind(or,lor,v,w),3)
lorbar <- sum(lor*w)/sum(w)
vlorbar <- 1/sum(w)

#summary odds ratio and 95% confidence interval
cbind(lorbar,vlorbar,ci(lorbar,vlorbar))
exp(c(lorbar,ci(lorbar,vlorbar)))

#SMR -- standard mortality ratio
E <- sum(p*R)/10^6
d <- sum(p*r)/10^6
SMR <- d/E
smr <- sum(p*r)/sum(p*R)
cbind(d,E,SMR,smr)
ci(log(smr),1/sum(d)+1/sum(D))
cbind(exp(ci(log(smr),1/sum(d)+1/sum(D))))

#smooth weighted average
x <- c(0.31,1.27,2.23,3.21,4.17,5.13,6.11,7.07,8.03,9.01)
y <- c(174.07,145.70,160.61,82.36,50.93,35.94,44.26,-21.78,86.27,242.46)

for(i in 1:length(x)) {
wt <- dnorm(x,x[i],1.0)
Y <- sum(wt*y)/sum(wt)
print(Y)
}

