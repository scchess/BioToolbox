#CHAPTER 19 -- Prediction performance

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}

#confidence intervals
ci(0.231,1/259+1/3595)
ci(0.109,0.031^2)
ci(0.0057,0.00213^2)
ci(.150,0.0416^2)

#data
sbp <- c(108,117,108,96,91,131,117,130,145,98,118,
111,86,94,125,91,133,105,133,110,99,118,
99,87,115,95,124,118,121,123,121,78,107,
119,89,100,158,117,106,101,105,118,118,104,
119,135,92,111,99,117)
dbp <- c(67,88,68,59,56,76,65,87,90,62,80,
72,52,60,70,57,66,63,76,67,51,79,
66,51,70,58,75,75,71,75,68,52,68,
74,58,70,79,66,72,56,67,73,70,71,
69,83,52,55,57,62)

#least squares estimated line
summary(lm(sbp~dbp))

#perpendicular
D <- .5*(var(dbp)-var(sbp))/cov(sbp,dbp)
B <- -D+sqrt(1+D^2)
B0 <- -D-sqrt(1+D^2)
A <- mean(sbp)-B*mean(dbp)
round(cbind(D,A,B,B0,-1/B0),3)
