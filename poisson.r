#CHAPTER 9 -- Poisson distribution

ci <- function(m,v) {
A <- m-1.960*sqrt(v)
B <- m+1.960*sqrt(v)
round(cbind(A,B),5)
}


#autism 
#data
cases <- c(202,255,333,442,571,716,803,996,1323)
births <- c(504853,534174,570976,613336,610701,602269,585761,568263,577113)
cbind(cases,births)

year <- 1987:1995
f <- glm(cases~year+offset(log(births)),family=poisson)
summary(f)
rates <- cases*10000/births
RATES <- fitted(f)*10000/births
CASES <- fitted(f)
X2 <- sum((cases-CASES)^2/CASES)
round(cbind(year,rates,RATES,cases,CASES),2)
pvalue <- 1-pchisq(X2,7)
round(cbind(X2,pvalue),3)
round(cbind(RATES[1],RATES[9],RATES[9]/RATES[1]),3)
round(cbind(rates[1],rates[9],rates[9]/rates[1]),3)

#analysis -- stomach cancer 45 to 85+
#data
deaths.m <- c(84,136,131,154,76)
pop.m <- c(2014272,1417097,775192,466018,161990)
base <-100000
rate.m <- base*deaths.m/pop.m
deaths.f <- c(52,69,77,121,92)
pop.f <- c(1967034,1474064,891864,642748,325940)
rate.f <- base*deaths.f/pop.f
round(cbind(deaths.m,pop.m,rate.m,deaths.f,pop.f,rate.f,rate.f/rate.m),3)
rate.ff <- sum(pop.f*rate.f)/sum(pop.f)
rate.mm <- sum(pop.f*rate.m)/sum(pop.f)
round(cbind(rate.ff,rate.mm,rate.ff/rate.mm),3)

#poisson model
sex <- sort(rep(0:1,5))
age <- c(45,55,65,75,85,45,55,65,75,85)
deaths <- c(deaths.m,deaths.f)
pop <- c(pop.m,pop.f)
cbind(deaths,pop,sex,age)
f <- glm(deaths~age+sex+offset(log(pop)),family=poisson)
summary(f)
round(rbind(deaths,fitted(f)),1)
#evaluation
X2 <- sum((deaths-fitted(f))^2/fitted(f))
pvalue <- 1-pchisq(X2,7)
cbind(X2,pvalue)

#ratios
#stomach cancer two young age groups
deaths <- c(44,162,44,135,22,44,12,58)
pop <- c(16569463,17170525,16652229,18283791,2829377,2644142,2981382,2975071)
base <- 1000000
rate <- base*deaths/pop
age <- c(0,1,0,1,0,1,0,1)
sex <- 1-c(1,1,0,0,1,1,0,0)
race <- c(0,0,0,0,1,1,1,1)
round(cbind(deaths,pop,rate,age,sex,race),2)

#stomach cancer ---age,sex and race
f <- glm(deaths~age+sex+race+offset(log(pop)),family=poisson)
summary(f)
DEATHS <- fitted(f)
RATES <- base*DEATHS/pop
round(rbind(DEATHS,RATES),1)
b <- summary(f)$coefficients[,1]
v <- (summary(f)$coefficients[,2])^2
round(c(1/122+1/399,1/272+1/249,1/136+1/385),7)
round(v[-1],7)

#log-variances, ratios and confidence intervals
cbind(round(v[-1],3),round(exp(b[-1]),3),round(exp(ci(b[-1],v[-1])),3))

#smoking pattern
#data
m <- matrix(c(98,1554,70,735,50,355,39,253),2,4)
m
smk <- c(0,10,20,30)
n <- m[1,]+m[2,]
p <- m[1,]/n
round(cbind(smk,p),3)
f <- glm(p~smk,weights=n,family=poisson)
summary(f)
round(rbind(log(fitted(f)),fitted(f),p),3)
#risk = 10 cigarettes/day
rr <- exp(10*coef(f)[2])
rr

#birth defects -- vitamin/ethnicity
x <- c(50,84,48,14,18,16,22,30,17,39,32,41)
d <- matrix(x,3,4)
d
chisq.test(d)
e <- outer(apply(d,1,sum),apply(d,2,sum))/sum(d)
round(e,1)
round(log(d),2)
round(log(e),2)

#detailed plot
matplot(1:3,log(d),type="l",xaxt="n",xlab="Categories (vitamin use)",
ylab="Log-frequencies",cex.lab=1.2)
matlines(1:3,log(e),lty=1)
title("Ethnicity by vitamin use")
text(1:3,c(2.6,2.6),"|")
text(c(2.5,2.5,2.5),c(4.4,3.8,3.3,3.0),c("white",
"Asian","Hispanic","African-American"))
text(c(1.15,2.15,2.85),c(2.6,2.6,2.6),c("always","during","never"))
legend(1,4.4,c("model","data"),lty=c(1,3),cex=0.75)

#analysis of incomplete table -- mother/daughter education
m <- matrix(c(84,142,43,100,106,48,67,77,92),3,3)
m
chisq.test(m)
rr <- factor(rep(1:3,3))
cc <- factor(sort(rr))
cbind(as.vector(m),rr,cc)
f <- glm(as.vector(m)~rr+cc,family=poisson)
summary(f)
#model estimated counts
round(matrix(fitted(f),3,3),3)
#data estimated counts
m0 <- outer(apply(m,1,sum),apply(m,2,sum))/sum(m)
round(m0,3)

#quasi independence
rr <- factor(c(2,3,1,3,1,2))
cc <- factor(c(1,1,2,2,3,3))
M <- c(142,43,100,48,67,77)
cbind(M,rr,cc)
f <- glm(M~rr+cc,family=poisson)
summary(f)
round(rbind(M,fitted(f)),3)
X2 <- sum((M-fitted(f))^2/fitted(f))
pvalue <- 1-pchisq(X2,1)
round(cbind(X2,pvalue),3)
