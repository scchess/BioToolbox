#CHAPTERS 23 -- Selection 

#mother/child pairs 
#data
A <- 9
B <- 14
C <- 16
D <- 25
E <- 18
F <- 15
G <- 23

#genetic analysis
o <- c(A,B,C,D,E,F,G)
n <- sum(o)
N  <- 3*n-D
p <- (3*A+2*B+2*C+D+E+F)/N
q <- 1-p
cbind(n,N,p)
e <- n*c(p^3,p^2*q,p^2*q,p*q,p*q^2,p*q^2,q^3)
round(rbind(o,e/n,e),3)
x2 <- sum((o-e)^2/e)
pvalue <- 1-pchisq(x2,6)
cbind(x2,pvalue)
v <- p*(1-p)/N
ci(p,v)

#Synder's ratios
#data
x11 <- 194
x12  <- 19
x21 <- 117
x22 <- 28
x31 <- 0
x32 <- 26
n <- sum(x11+x12+x21+x22+x31+x32)
q <- sqrt((x12+x22+x32)/n)
p <- 1-q
cbind(q,p)

#calculation of Synder's ratios
s2 <- x12/(x11+x12)
s1 <- x22/(x21+x22)
c(s1,s2)
s2/s1^2

#genetic chi-square analysis
e11 <- p^2*(1+2*q)
e12 <- p^2*q^2
e21 <- 2*p*q^2
e22 <- 2*p*q^3
E11 <- (x11+x12)*e11/(e11+e12)
E12 <- (x11+x12)*e12/(e11+e12)
E21 <- (x21+x22)*e21/(e21+e22)
E22 <- (x21+x22)*e22/(e21+e22)

o <- c(x11,x12,x21,x22)
e <- c(E11,E12,E21,E22)
cbind(round(e,1),o)
xx <- sum((o-e)^2/e)
pvalue <- 1-pchisq(xx,2)
cbind(xx,pvalue)


#Selection
#balanced polymorphism
s1 <- 0.2
s2 <- 0
s3 <- 0.6
p <- 0.2
#s=1.0 -- aa fails to survive 
s1 <- 0
s2 <- 0
s3 <- 1
p <- 0.3
#s=0.5 or 0.1 -- survival probability = 0.5 or 0.1
s1 <- 0
s2 <- 0
s3 <- 0.5
s3 <- 0.1
p <- 0.3

#analysis
q <- 1-p
for(ii in 1:40) {
D <- p^2*(1-s1)
H <- 2*p*q*(1-s2)
R <- q^2*(1-s3)
p <- (D+.5*H)/(D+H+R)
q <- 1-p
if(ii==1|ii==2|ii==3|ii==10|ii==20|ii==40) print(round(cbind(ii,D,H,R,p),3))
}

#Mr.A versus mr.B
a <- 100
b <- 10
pa <- 0.5
pb <- 0.1

for(i in 1:15) {
a0 <- a-a*pa+b*pb
b0 <- b+a*pa-b*pb
a.to.b <- a*pa
b.to.a <- b*pb
ratio <- a/b
print(round(cbind(i,a0,b0,a.to.b,b.to.a,ratio),2))
a <- a0
b <- b0
}
