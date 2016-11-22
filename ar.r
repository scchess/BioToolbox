#CHAPTER 20 -- Attributable risk

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#usual estimation form a 2 by 2 table
a <- 159
b <- 1343
c <- 98
d <- 1553
n <- a+b+c+d
P <- (a+c)/n
p <- c/(c+d)
rr <- (a/(a+b))/(c/(c+d))
e <- (a+b)/n
round(cbind(P,p,rr,e),3)

#attributable risk estimates
#version 1
(P-p)/P
#version 2
(a*d-b*c)/((a+c)*(c+d))
#version 3
e*(rr-1)/(e*(rr-1)+1)
#version 4
1-1/(e*rr+(1-e)*1.0)
#version 5
(a/(a+c))*((rr-1)/rr)

#confidence intervals -- log(1-ar) and ar
ar <- (P-p)/P
V <- (b+ar*(a+d))/(n*c)
ci(log(1-ar),V)
rev(1-exp(ci(log(1-ar),V)))

#age-strata summary
d0 <- c(5,15,21,15)
P0 <- c(100,200,350,600)
d1 <- c(10,36,54,40)
P1 <- c(100,200,300,500)
D <- d1/sum(d1)
rr <- (d1/P1)/(d0/P0)
e <- sum(P1+d1)/(sum(P0+d0)+sum(P1+d1))
ar <- e*(rr-1)/(e*(rr-1)+1)
round(cbind(D,rr,ar),3)
arbar <- sum(D*ar)
arbar


#chd and smoking
chd <- c(59,37,100,61)
nochd <- c(271,325,1072,1228)
r <- (chd/(chd+nochd))/(chd[4]/(chd[4]+nochd[4]))
e <- (chd+nochd)/sum(chd+nochd)
round(cbind(e,r),3)
P <- sum(chd)/sum(chd+nochd)
p <- (chd[4]/(chd[4]+nochd[4]))
ar1 <- (P-p)/P
ar2 <- 1-1/sum(e*r)
round(cbind(P,p,ar1,ar2),3)

#data -- chd and smoking at four level of cholesterol
a <- c(14,17,47,81)
b <- c(333,248,381,381)
c <- c(17,12,27,41)
d <- c(467,343,374,369)
n <- a+b+c+d
#estimates
D <- (a+c)/(sum(a+c))
p <- c/(c+d)
P <- (a+c)/n
ar <- (P-p)/P
v <- ((1-ar)^2*(b+ar*(a+d)))/(n*c)
round(cbind(D,p,P,ar,v),3)

#weighted average -- estimated attributable risk  
#version 1
w <- 1/v
V <- 1/sum(w)
arbar <- sum(w*ar)/sum(w)
arbar
ci(arbar,V)
#version 2
arbar0 <- sum(D*ar)
arbar0
ci(arbar0,V)

#simple misclassification model
A <- c(1.0,0.95,0.90,0.80)
a <- 159
b <- 1343
c <- 98
d <- 1553
n <- a+b+c+d
p <- (c*A)/(c*A+d)
P <- (a+c)/n
ar <- (P-p)/P
round(ar,3)
