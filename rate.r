#CHAPTER 16 -- Rates 

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
cbind(A,B)
}

#confidence intervals
d <- 234
pop <- 11139194
r <- d/pop
ci(log(r),1/d)
exp(ci(log(r),1/d))*100000

d0 <- 549
pop0 <- 63818574
r0 <- d0/pop0
rr <- r/r0
cbind(10^5*r,10^5*r0,rr)
v <- 1/234+1/549
ci(log(rr),v)
exp(ci(log(rr),v))

#Mens Health Study data
#data: non-smokers
smk0 <- c(2,42,27,22,26,16,31,37,15,30,12,5,80,29,13,1,14)
c.smk0 <- c(0,0,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1)

#data: non-smokers
d0 <- sum(c.smk0)
r0 <- d0/sum(smk0) 
v0 <- 1/d0
cbind(r0,log(r0),v0)
#ci --- log-rate
ci(log(r0),v0)
#ci --- rate
exp(ci(log(r0),v0))
rev(1/exp(ci(log(r0),v0)))
#mean survival time
1/r0
median
m0 <- log(2)/r0
m0
#ci --- log-median
ci(log(m0),v0)
#ci --- median
exp(ci(log(m0),v0))
rev(log(2)/exp(ci(log(r0),v0)))

#data: smokers
smk <- c(21,13,17,8,23,18)
c.smk <- c(0,1,1,1,1,1)

#smokers
d <- sum(c.smk)
r <- d/sum(smk) 
v <- 1/d
cbind(r,log(r),v)
#mean survival time
1/r
#median
m <- log(2)/r
m
#ci --- log-rate
ci(log(r),v)
#ci --- rate
exp(ci(log(r),v))
rev(1/exp(ci(log(r),v)))
#ci --- log-median
ci(log(m),v)
#ci --- median
exp(ci(log(m),v))
rev(log(2)/exp(ci(log(r),v)))
