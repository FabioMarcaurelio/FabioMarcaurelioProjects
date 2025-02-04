rm(list=ls())


library(sandwich)
library(bucky)
library(lmtest)

### ADVANCED ECONOMETRICS ###
### HOMEWORK 2 ###
### Fabio Marcaurelio (599027)###


################################
### FOR COMMENTS SEE PDF FILE###
################################


### simulate R=1000 A n=50 ###

set.seed(599027)
beta1=1
beta2=1
R=1000
coefA1_50<-numeric(R)
coefA2_50<-numeric(R)
varbetaA1H_50<-numeric(R)
varbetaA2H_50<-numeric(R)
varbetaA1R_50<-numeric(R)
varbetaA2R_50<-numeric(R)
sqA1<-numeric(R)
sqA2<-numeric(R)
lowerA1_50<- numeric(R)
upperA1_50<- numeric(R)
lowerA2_50<- numeric(R)
upperA2_50<-numeric(R)
lowerA1_50R<-numeric(R)
upperA1_50R<-numeric(R)
lowerA2_50R<-numeric(R)
upperA2_50R<-numeric(R)

for (i in 1:R) {
  nA1=50
  uA1<-rnorm(nA1,0,1)
  x1<-rnorm(nA1,0,1)
  x2<-rnorm(nA1,0,1)
  x3<-rnorm(nA1,0,1)
  xA1<-x1^2+x2^2+x3^2
  
  yA1<-beta1+beta2*xA1+uA1
  olsA1<-lm(yA1~xA1)
  sampleA1<-olsA1

  coefA1_50[i]<-olsA1$coefficients[1]
  coefA2_50[i]<-olsA1$coefficients[2]
  
  homoA1<-vcov(sampleA1)
  homoA1
  varbetaA1H_50[i]<-vcov(sampleA1)[1,1]
  varbetaA2H_50[i]<-vcov(sampleA1)[2,2]
  varbetaA1R_50[i]<-vcovHC(sampleA1,type="HC0")[1,1]
  varbetaA2R_50[i]<-vcovHC(sampleA1,type="HC0")[2,2]
  
  sqA1<-sqrt(vcov(sampleA1)[1,1])
  sqA2<-sqrt(vcov(sampleA1)[2,2])
  vcmrA1 <- vcovHC(sampleA1, type="HC0")
  vcmrA1
  
  lowerA1_50[i] <- (coefA1_50[i]) - qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA1H_50[i])
  upperA1_50[i] <- (coefA1_50[i]) + qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA1H_50[i])
  
  lowerA2_50[i] <- (coefA2_50[i]) - qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA2H_50[i])
  upperA2_50[i] <- (coefA2_50[i]) + qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA2H_50[i])
  
  lowerA1_50R[i] <- (coefA1_50[i]) - qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA1R_50[i])
  upperA1_50R[i] <- (coefA1_50[i]) + qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA1R_50[i])
  
  lowerA2_50R[i] <- (coefA2_50[i]) - qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA2R_50[i])
  upperA2_50R[i] <- (coefA2_50[i]) + qt(0.975, df = sampleA1$df[1]) * sqrt(varbetaA2R_50[i])
}

### 0.95 coverage test for t under homoskedasticity and robust to heteroskedasticity  ###

### 0.95 covarage for beta=1 homoskedasticity ###

CIA1_50 <- cbind(lowerA1_50, upperA1_50)
mean(CIA1_50[, 1] <= 1 & 1 <= CIA1_50[, 2])
CIA2_50 <- cbind(lowerA2_50, upperA2_50)
mean(CIA2_50[, 1] <= 1 & 1 <= CIA2_50[, 2])

### 0.95 covarage for beta=0 homoskedasticity ###

CIA1_50zero <- cbind(lowerA1_50, upperA1_50)
mean(CIA1_50zero[, 1] <= 0 & 0 <= CIA1_50[, 2])
CIA2_50zero <- cbind(lowerA2_50, upperA2_50)
mean(CIA2_50zero[, 1] <= 0 & 0 <= CIA2_50zero[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty ###

CIA1_50R <- cbind(lowerA1_50R, upperA1_50R)
mean(CIA1_50R[, 1] <= 1 & 1 <= CIA1_50[, 2])
CIA2_50R <- cbind(lowerA2_50R, upperA2_50R)
mean(CIA2_50R[, 1] <= 1 & 1 <= CIA2_50R[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty ###

CIA1_50zeroR <- cbind(lowerA1_50R, upperA1_50R)
mean(CIA1_50zeroR[, 1] <= 0 & 0 <= CIA1_50R[, 2])
CIA2_50zeroR <- cbind(lowerA2_50R, upperA2_50R)
mean(CIA2_50zeroR[, 1] <= 0 & 0 <= CIA2_50zeroR[, 2])

### plot Y ###

plot(xA1, yA1,
     ylab="yA1", xlab="xA1")
abline(sampleA1, col="blue")

cov(xA1,uA1)

### comapre variancee-covariance matrix under hoskedastic and robust to heteroskedastic ###

varA1H_50<-mean(varbetaA1H_50)
varA2H_50<-mean(varbetaA2H_50)
varA1R_50<-mean(varbetaA1R_50)
varA2R_50<-mean(varbetaA2R_50)
c(varA1H_50,varA1R_50)
c(varA2H_50,varA2R_50)

homoA1
vcmrA1
sqA1
sqA2

### t-test under beta=1 and beta=0 ###

summary(sampleA1)
robustolsA1<-robust.summary(sampleA1, vcov = vcovHC(sampleA1), type = "HC3")

(tstatA1 <- (robustolsA1$coefficients[2,1]-1)/robustolsA1$coefficients[2,2])
(pvalueA1 <- (pt(tstatA1, df = 50-2, lower.tail = TRUE)*2))
(pvalueNA1 <- (pnorm(tstatA1, lower.tail = FALSE)*2))

(tstatA10 <- (robustolsA1$coefficients[2,1])/robustolsA1$coefficients[2,2])
(pvalueA10 <- (pt(tstatA10, df = 50-2, lower.tail = TRUE)*2))
(pvalueA10 <- (pnorm(tstatA10, lower.tail = FALSE)*2))


### simulate R=1000 A n=500 ###

set.seed(599027)
R=1000
coefA1_500<-numeric(R)
coefA2_500<-numeric(R)
varbetaA1R_500<-numeric(R)
varbetaA2R_500<-numeric(R)
varbetaA1H_500<-numeric(R)
varbetaA2H_500<-numeric(R)
lowerA1_500<- numeric(R)
upperA1_500<- numeric(R)
lowerA2_500<- numeric(R)
upperA2_500<-numeric(R)
lowerA1_500R<-numeric(R)
upperA1_500R<-numeric(R)
lowerA2_500R<-numeric(R)
upperA2_500R<-numeric(R)

for (i in 1:R) {
  nA2=500
  uA2<-rnorm(nA2,0,1)
  x1<-rnorm(nA2,0,1)
  x2<-rnorm(nA2,0,1)
  x3<-rnorm(nA2,0,1)
  xA2<-x1^2+x2^2+x3^2
  
  yA2<-beta1+beta2*xA2+uA2
  olsA2<-lm(yA2~xA2)
  sampleA2<-olsA2
  
  coefA1_500[i]<-olsA2$coefficients[1]
  coefA2_500[i]<-olsA2$coefficients[2]
  
  homoA2<-vcov(sampleA2)
  homoA2
  
  varbetaA1H_500[i]<-vcov(sampleA2)[1,1]
  varbetaA2H_500[i]<-vcov(sampleA2)[2,2]
  varbetaA1R_500[i]<-vcovHC(sampleA2,type="HC0")[1,1]
  varbetaA2R_500[i]<-vcovHC(sampleA2,type="HC0")[2,2]
  
  sqA1_500<-sqrt(vcov(sampleA2)[1,1])
  sqA2_500<-sqrt(vcov(sampleA2)[2,2])
  vcmrA2 <- vcovHC(sampleA2, type="HC0")
  vcmrA2
  robust.summary(sampleA2, vcov = vcovHC(sampleA2), type = "HC3")
  sqrt(diag(vcovHC(sampleA2, type="HC3")))
  
  lowerA1_500[i] <- (coefA1_500[i]) - 1.96 * sqrt(varbetaA1H_500[i])
  upperA1_500[i] <- (coefA1_500[i]) + 1.96 * sqrt(varbetaA1H_500[i]) 
  
  lowerA2_500[i] <- (coefA2_500[i]) - 1.96 * sqrt(varbetaA2H_500[i])
  upperA2_500[i] <- (coefA2_500[i]) + 1.96 * sqrt(varbetaA2H_500[i]) 
  
  lowerA1_500R[i] <- (coefA1_500[i]) - 1.96 * sqrt(varbetaA1R_500[i])
  upperA1_500R[i] <- (coefA1_500[i]) + 1.96 * sqrt(varbetaA1R_500[i]) 
  
  lowerA2_500R[i] <- (coefA2_500[i]) - 1.96 * sqrt(varbetaA2R_500[i])
  upperA2_500R[i] <- (coefA2_500[i]) + 1.96 * sqrt(varbetaA2R_500[i]) 
  
}

### 0.95 coverage test for t under homoskedasticity and robust to heteroskedasticity  ###

### 0.95 covarage for beta=1 homoskedasticity ###

CIA1_500 <- cbind(lowerA1_500, upperA1_500)
mean(CIA1_500[, 1] <= 1 & 1 <= CIA1_500[, 2])
CIA2_500 <- cbind(lowerA2_500, upperA2_500)
mean(CIA2_500[, 1] <= 1 & 1 <= CIA2_500[, 2])

### 0.95 covarage for beta=0 homoskedasticity ###

CIA1_500zero <- cbind(lowerA1_500, upperA1_500)
mean(CIA1_500zero[, 1] <= 0 & 0 <= CIA1_500zero[, 2])
CIA2_500zero <- cbind(lowerA2_500, upperA2_500)
mean(CIA2_500zero[, 1] <= 0 & 0 <= CIA2_500zero[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty ###

CIA1_500R <- cbind(lowerA1_500R, upperA1_500R)
mean(CIA1_500R[, 1] <= 1 & 1 <= CIA1_500R[, 2])
CIA2_500R <- cbind(lowerA2_500R, upperA2_500R)
mean(CIA2_500R[, 1] <= 1 & 1 <= CIA2_500R[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty ###

CIA1_500zeroR <- cbind(lowerA1_500R, upperA1_500R)
mean(CIA1_500zeroR[, 1] <= 0 & 0 <= CIA1_500zeroR[, 2])
CIA2_500zeroR <- cbind(lowerA2_500R, upperA2_500R)
mean(CIA2_500zeroR[, 1] <= 0 & 0 <= CIA2_500zeroR[, 2])


### plot Y ###

plot(xA2, yA2,
     ylab="yA2", xlab="xA2")
abline(sampleA2, col="blue")

cov(xA2,uA2)

### comapre variancee-covariance matrix under hoskedastic and robust to heteroskedastic ###

varA1H_500<-mean(varbetaA1H_500)
varA2H_500<-mean(varbetaA2H_500)
varA1R_500<-mean(varbetaA1R_500)
varA2R_500<-mean(varbetaA2R_500)
c(varA1H_500,varA1R_500)
c(varA2H_500,varA2R_500)

homoA2
vcmrA2
sqA1_500
sqA2_500

### t-test under beta=1 and beta=0 ###

summary(sampleA2)
robustolsA2<-robust.summary(sampleA2, vcov = vcovHC(sampleA2), type = "HC3")

(tstatA2 <- (robustolsA2$coefficients[2,1]-1)/robustolsA2$coefficients[2,2])
(pvalueA2 <- (pt(tstatA2, df = 50-2, lower.tail = TRUE)*2))
(pvalueNA2 <- (pnorm(tstatA2, lower.tail = FALSE)*2))

(tstatA20 <- (robustolsA2$coefficients[2,1])/robustolsA2$coefficients[2,2])
(pvalueA20 <- (pt(tstatA20, df = 50-2, lower.tail = TRUE)*2))
(pvalueA20 <- (pnorm(tstatA20, lower.tail = FALSE)*2))

#################################################################################################

################################
### FOR COMMENTS SEE PDF FILE###
################################

### simulate R=1000 B n=50###

set.seed(599027)
R=1000
coefB1_50<-numeric(R)
coefB2_50<-numeric(R)
varbetaB1R_50<-numeric(R)
varbetaB2R_50<-numeric(R)
varbetaB1H_50<-numeric(R)
varbetaB2H_50<-numeric(R)
lowerB1_50<- numeric(R)
upperB1_50<- numeric(R)
lowerB2_50<- numeric(R)
upperB2_50<-numeric(R)
lowerB1_50R<-numeric(R)
upperB1_50R<-numeric(R)
lowerB2_50R<-numeric(R)
upperB2_50R<-numeric(R)
lowerB1_50A<-numeric(R)
upperB1_50A<-numeric(R)
lowerB2_50A<-numeric(R)
upperB2_50A<-numeric(R)
lowerB1R_50A<-numeric(R)
upperB1R_50A<-numeric(R)
lowerB2R_50A<-numeric(R)
upperB2R_50A<-numeric(R)
for (i in 1:R) {
  nB1=50
  uB1<-runif(nB1,-3,3)
  x1<-rnorm(nB1,0,1)
  x2<-rnorm(nB1,0,1)
  x3<-rnorm(nB1,0,1)
  xB1<-x1^2+x2^2+x3^2
  
  yB1<-beta1+beta2*xB1+uB1
  olsB1<-lm(yB1~xB1)
  sampleB1<-olsB1
  coefB1_50[i]<-sampleB1$coefficients[1]
  coefB2_50[i]<-sampleB1$coefficients[2]
  
  homoB1<-vcov(sampleB1)
  homoB1
  
  varbetaB1H_50[i]<-vcov(sampleB1)[1,1]
  varbetaB2H_50[i]<-vcov(sampleB1)[2,2]
  varbetaB1R_50[i]<-vcovHC(sampleB1,type="HC0")[1,1]
  varbetaB2R_50[i]<-vcovHC(sampleB1,type="HC0")[2,2]
  
  sqB1_50<-sqrt(vcov(sampleB1)[1,1])
  aqB1_50<-sqrt(vcov(sampleB1)[2,2])
  vcmrB_50 <- vcovHC(sampleB1, type="HC0")
  vcmrB_50
  
  lowerB1_50A[i] <- (coefB1_50[i]) - 1.96 * sqrt(varbetaB1H_50[i])
  upperB1_50A[i] <- (coefB1_50[i]) + 1.96 * sqrt(varbetaB1H_50[i]) 
  
  lowerB2_50A[i] <- (coefB2_50[i]) - 1.96 * sqrt(varbetaB2H_50[i])
  upperB2_50A[i] <- (coefB2_50[i]) + 1.96 * sqrt(varbetaB2H_50[i]) 
  
  lowerB1R_50A[i] <- (coefB1_50[i]) - 1.96 * sqrt(varbetaB1R_50[i])
  upperB1R_50A[i] <- (coefB1_50[i]) + 1.96 * sqrt(varbetaB1R_50[i]) 
  
  lowerB2R_50A[i] <- (coefB2_50[i]) - 1.96 * sqrt(varbetaB2R_50[i])
  upperB2R_50A[i] <- (coefB2_50[i]) + 1.96 * sqrt(varbetaB2R_50[i])
  
  lowerB1_50[i] <- (coefB1_50[i]) - qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB1H_50[i])
  upperB1_50[i] <- (coefB1_50[i]) + qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB1H_50[i])
  
  lowerB2_50[i] <- (coefB2_50[i]) - qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB2H_50[i])
  upperB2_50[i] <- (coefB2_50[i]) + qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB2H_50[i])
  
  lowerB1_50R[i] <- (coefB1_50[i]) - qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB1R_50[i])
  upperB1_50R[i] <- (coefB1_50[i]) + qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB1R_50[i])
  
  lowerB2_50R[i] <- (coefB2_50[i]) - qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB2R_50[i])
  upperB2_50R[i] <- (coefB2_50[i]) + qt(0.975, df = sampleB1$df[1]) * sqrt(varbetaB2R_50[i])
}

### 0.95 coverage test for z and t under homoskedasticity and robust to heteroskedasticity ###

### 0.95 covarage for beta=1 homoskedasticity T-STUDENT ###

CIB1_50 <- cbind(lowerB1_50, upperB1_50)
mean(CIB1_50[, 1] <= 1 & 1 <= CIB1_50[, 2])
CIB2_50 <- cbind(lowerB2_50, upperB2_50)
mean(CIB2_50[, 1] <= 1 & 1 <= CIB2_50[, 2])

### 0.95 covarage for beta=0 homoskedasticity T-STUDENT ###

CIB1_50zero <- cbind(lowerB1_50, upperB1_50)
mean(CIB1_50zero[, 1] <= 0 & 0 <= CIB1_50zero[, 2])
CIB2_50zero <- cbind(lowerB2_50, upperB2_50)
mean(CIB2_50zero[, 1] <= 0 & 0 <= CIB2_50zero[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty T-STUDENT ###

CIB1_50R <- cbind(lowerB1_50R, upperB1_50R)
mean(CIB1_50R[, 1] <= 1 & 1 <= CIB1_50R[, 2])
CIB2_50R <- cbind(lowerB2_50R, upperB2_50R)
mean(CIB2_50R[, 1] <= 1 & 1 <= CIB2_50R[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty T-STUDENT ###

CIB1_50zeroR <- cbind(lowerB1_50R, upperB1_50R)
mean(CIB1_50zeroR[, 1] <= 0 & 0 <= CIB1_50zeroR[, 2])
CIB2_50zeroR <- cbind(lowerB2_50R, upperB2_50R)
mean(CIB2_50zeroR[, 1] <= 0 & 0 <= CIB2_50zeroR[, 2])

### 0.95 covarage for beta=1 homoskedasticity NORMAL ###

CIB1_50A <- cbind(lowerB1_50A, upperB1_50A)
mean(CIB1_50A[, 1] <=  1 &  1 <= CIB1_50A[, 2])
CIB2_50A <- cbind(lowerB2_50A, upperB2_50A)
mean(CIB2_50A[, 1] <= 1 & 1 <= CIB2_50A[, 2])

### 0.95 covarage for beta=0 homoskedasticity NORMAL ###

CIB1_50zeroA <- cbind(lowerB1_50A, upperB1_50A)
mean(CIB1_50zeroA[, 1] <=  0 &  0 <= CIB1_50zeroA[, 2])
CIB2_50zeroA <- cbind(lowerB2_50A, upperB2_50A)
mean(CIB2_50zeroA[, 1] <= 0 & 0 <= CIB2_50zeroA[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty NORMAL ###

CIB1R_50A <- cbind(lowerB1R_50A, upperB1R_50A)
mean(CIB1R_50A[, 1] <=  1 &  1 <= CIB1R_50A[, 2])
CIB2R_50A <- cbind(lowerB2R_50A, upperB2R_50A)
mean(CIB2R_50A[, 1] <= 1 & 1 <= CIB2R_50A[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty NORMAL ###

CIB1R_50zeroA <- cbind(lowerB1R_50A, upperB1R_50A)
mean(CIB1R_50zeroA[, 1] <=  0 &  0 <= CIB1R_50zeroA[, 2])
CIB2R_50zeroA <- cbind(lowerB2R_50A, upperB2R_50A)
mean(CIB2R_50zeroA[, 1] <= 0 & 0 <= CIB2R_50zeroA[, 2])

### plot Y ###

plot(xB1, yB1,
     ylab="yB1", xlab="xB1")
abline(sampleB1, col="blue")

cov(xB1,uB1)

### comapre variancee-covariance matrix under hoskedastic and robust to heteroskedastic ###

varB1H_50<-mean(varbetaB1H_50)
varB2H_50<-mean(varbetaB2H_50)
varB1R_50<-mean(varbetaB1R_50)
varB2R_50<-mean(varbetaB2R_50)
c(varB1H_50,varB1R_50)
c(varB2H_50,varB2R_50)

homoB1
vcmrB_50
sqA1_50
sqA2_50

### t-test under beta=1 and beta=0 ###
 
summary(sampleB1)
robustolsB1<-robust.summary(sampleB1, vcov = vcovHC(sampleB1), type = "HC3")

(tstatB1 <- (robustolsB1$coefficients[2,1]-1)/robustolsB1$coefficients[2,2])
(pvalueB1 <- (pt(tstatB1, df = 50-2, lower.tail = TRUE)*2))
(pvalueNB1 <- (pnorm(tstatB1, lower.tail = FALSE)*2))

(tstatB10 <- (robustolsB1$coefficients[2,1])/robustolsB1$coefficients[2,2])
(pvalueB10 <- (pt(tstatB10, df = 50-2, lower.tail = TRUE)*2))
(pvalueB10 <- (pnorm(tstatB10, lower.tail = FALSE)*2))


### simulate R=1000 B n=500###

set.seed(599027)
R=1000
coefB1_500<-numeric(R)
coefB2_500<-numeric(R)
varbetaB1R_500<-numeric(R)
varbetaB2R_500<-numeric(R)
varbetaB1H_500<-numeric(R)
varbetaB2H_500<-numeric(R)
lowerB1_500<- numeric(R)
upperB1_500<- numeric(R)
lowerB2_500<- numeric(R)
upperB2_500<-numeric(R)
lowerB1R_500<-numeric(R)
upperB1R_500<-numeric(R)
lowerB2R_500<-numeric(R)
upperB2R_500<-numeric(R)
lowerB1_500A<-numeric(R)
upperB1_500A<-numeric(R)
lowerB2_500A<-numeric(R)
upperB2_500A<-numeric(R)
lowerB1R_500A<-numeric(R)
upperB1R_500A<-numeric(R)
lowerB2R_500A<-numeric(R)
upperB2R_500A<-numeric(R)
for (i in 1:R) {
  nB2=500
  uB2<-runif(nB2,-3,3)
  x1<-rnorm(nB2,0,1)
  x2<-rnorm(nB2,0,1)
  x3<-rnorm(nB2,0,1)
  xB2<-x1^2+x2^2+x3^2
  
  yB2<-beta1+beta2*xB2+uB2
  olsB2<-lm(yB2~xB2)
  sampleB2<-olsB2
  coefB1_500[i]<-olsB2$coefficients[1]
  coefB2_500[i]<-olsB2$coefficients[2]
  
  homoB2<-vcov(sampleB2)
  homoB2
  
  varbetaB1H_500[i]<-vcov(sampleB2)[1,1]
  varbetaB2H_500[i]<-vcov(sampleB2)[2,2]
  varbetaB1R_500[i]<-vcovHC(sampleB2,type="HC0")[1,1]
  varbetaB2R_500[i]<-vcovHC(sampleB2,type="HC0")[2,2]
  
  sqB1_500<-sqrt(vcov(sampleB2)[1,1])
  sqB2_500<-sqrt(vcov(sampleB2)[2,2])
  vcmrB_500<- vcovHC(sampleB2, type="HC0")
  vcmrB_500
  
  lowerB1_500[i] <- (coefB1_500[i]) - 1.96 * sqrt(varbetaB1H_500[i])
  upperB1_500[i] <- (coefB1_500[i]) + 1.96 * sqrt(varbetaB1H_500[i]) 
  
  lowerB2_500[i] <- (coefB2_500[i]) - 1.96 * sqrt(varbetaB2H_500[i])
  upperB2_500[i] <- (coefB2_500[i]) + 1.96 * sqrt(varbetaB2H_500[i]) 
  
  lowerB1R_500[i] <- (coefB1_500[i]) - 1.96 * sqrt(varbetaB1R_500[i])
  upperB1R_500[i] <- (coefB1_500[i]) + 1.96 * sqrt(varbetaB1R_500[i]) 
  
  lowerB2R_500[i] <- (coefB2_500[i]) - 1.96 * sqrt(varbetaB2R_500[i])
  upperB2R_500[i] <- (coefB2_500[i]) + 1.96 * sqrt(varbetaB2R_500[i]) 
  
  lowerB1_500A[i] <- (coefB1_500[i]) - qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB1H_500[i])
  upperB1_500A[i] <- (coefB1_500[i]) + qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB1H_500[i]) 
  
  lowerB2_500A[i] <- (coefB2_500[i]) - qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB2H_500[i])
  upperB2_500A[i] <- (coefB2_500[i]) + qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB2H_500[i]) 
  
  lowerB1R_500A[i] <- (coefB1_500[i]) - qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB1R_500[i])
  upperB1R_500A[i] <- (coefB1_500[i]) + qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB1R_500[i]) 
  
  lowerB2R_500A[i] <- (coefB2_500[i]) - qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB2R_500[i])
  upperB2R_500A[i] <- (coefB2_500[i]) + qt(0.975, df = sampleB2$df[1]) * sqrt(varbetaB2R_500[i]) 
  
}

### 0.95 coverage test for z and t under homoskedasticity and robust to heteroskedasticity ###

### 0.95 covarage for beta=1 homoskedasticity NORMAL ###

CIB1_500 <- cbind(lowerB1_500, upperB1_500)
mean(CIB1_500[, 1] <=  1 &  1 <= CIB1_500[, 2])
CIB2_500 <- cbind(lowerB2_500, upperB2_500)
mean(CIB2_500[, 1] <= 1 & 1 <= CIB2_500[, 2])

### 0.95 covarage for beta=0 homoskedasticity NORMAL ###

CIB1_500zero <- cbind(lowerB1_500, upperB1_500)
mean(CIB1_500zero[, 1] <=  0 &  0 <= CIB1_500zero[, 2])
CIB2_500zero <- cbind(lowerB2_500, upperB2_500)
mean(CIB2_500zero[, 1] <= 0 & 0 <= CIB2_500zero[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty NORMAL ###

CIB1R_500 <- cbind(lowerB1R_500, upperB1R_500)
mean(CIB1R_500[, 1] <=  1 &  1 <= CIB1R_500[, 2])
CIB2R_500 <- cbind(lowerB2R_500, upperB2R_500)
mean(CIB2R_500[, 1] <= 1 & 1 <= CIB2R_500[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty NORMAL ###

CIB1R_500zero <- cbind(lowerB1R_500, upperB1R_500)
mean(CIB1R_500zero[, 1] <=  0 &  0 <= CIB1R_500zero[, 2])
CIB2R_500zero <- cbind(lowerB2R_500, upperB2R_500)
mean(CIB2R_500zero[, 1] <= 0 & 0 <= CIB2R_500zero[, 2])

### 0.95 covarage for beta=1 homoskedasticity T-STUDENT ###

CIB1_500A <- cbind(lowerB1_500A, upperB1_500A)
mean(CIB1_500A[, 1] <=  1 &  1 <= CIB1_500A[, 2])
CIB2_500A <- cbind(lowerB2_500A, upperB2_500A)
mean(CIB2_500A[, 1] <= 1 & 1 <= CIB2_500A[, 2])

### 0.95 covarage for beta=0 homoskedasticity T-STUDENT ###

CIB1_500zeroA <- cbind(lowerB1_500A, upperB1_500A)
mean(CIB1_500zeroA[, 1] <=  0 &  0 <= CIB1_500zeroA[, 2])
CIB2_500zeroA <- cbind(lowerB2_500A, upperB2_500A)
mean(CIB2_500zeroA[, 1] <= 0 & 0 <= CIB2_500zeroA[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty T-STUDENT ###

CIB1R_500A <- cbind(lowerB1R_500A, upperB1R_500A)
mean(CIB1R_500A[, 1] <=  1 &  1 <= CIB1R_500A[, 2])
CIB2R_500A <- cbind(lowerB2R_500A, upperB2R_500A)
mean(CIB2R_500A[, 1] <= 1 & 1 <= CIB2R_500A[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty T-STUDENT ###

CIB1R_500zeroA <- cbind(lowerB1R_500A, upperB1R_500A)
mean(CIB1R_500zeroA[, 1] <=  0 &  0 <= CIB1R_500zeroA[, 2])
CIB2R_500zeroA <- cbind(lowerB2R_500A, upperB2R_500A)
mean(CIB2R_500zeroA[, 1] <= 0 & 0 <= CIB2R_500zeroA[, 2])

### plot Y ###

plot(xB2, yB2,
     ylab="yB2", xlab="xB2")
abline(sampleB2, col="blue")

cov(xB2,uB2)

### comapre variancee-covariance matrix under hoskedastic and robust to heteroskedastic ###

varB1H_500<-mean(varbetaB1H_500)
varB2H_500<-mean(varbetaB2H_500)
varB1R_500<-mean(varbetaB1R_500)
varB2R_500<-mean(varbetaB2R_500)
c(varB1H_500,varB1R_500)
c(varB2H_500,varB2R_500)

homoB2
vcmrB_500

sqB1_500
sqB2_500

### t-test under beta=1 and beta=0 ###

summary(sampleB2)
robustolsB2<-robust.summary(sampleB2, vcov = vcovHC(sampleB2), type = "HC3")

(tstatB2 <- (robustolsB2$coefficients[2,1]-1)/robustolsB2$coefficients[2,2])
(pvalueB2 <- (pt(tstatB2, df = 500-2, lower.tail = TRUE)*2))
(pvalueNB2 <- (pnorm(tstatB2, lower.tail = FALSE)*2))

(tstatB20 <- (robustolsB2$coefficients[2,1])/robustolsB2$coefficients[2,2])
(pvalueB20 <- (pt(tstatB20, df = 500-2, lower.tail = TRUE)*2))
(pvalueB20 <- (pnorm(tstatB20, lower.tail = FALSE)*2))

#################################################################################################

################################
### FOR COMMENTS SEE PDF FILE###
################################

### simulate R=1000 C n=100 ###


set.seed(599027)
R=1000
coefC1_100<-numeric(R)
coefC2_100<-numeric(R)
varbetaC1R_100<-numeric(R)
varbetaC2R_100<-numeric(R)
varbetaC1H_100<-numeric(R)
varbetaC2H_100<-numeric(R)
lowerC1_100<- numeric(R)
upperC1_100<- numeric(R)
lowerC2_100<- numeric(R)
upperC2_100<-numeric(R)
lowerC1R_100<- numeric(R)
upperC1R_100<- numeric(R)
lowerC2R_100<- numeric(R)
upperC2R_100<-numeric(R)
lowerC1_100A<-numeric(R)
upperC1_100A<-numeric(R)
lowerC2_100A<-numeric(R)
upperC2_100A<-numeric(R)
lowerC1R_100A<-numeric(R)
upperC1R_100A<-numeric(R)
lowerC2R_100A<-numeric(R)
upperC2R_100A<-numeric(R)
for (i in 1:R) {
  nC=100
  x1<-rnorm(nC,0,1)
  x2<-rnorm(nC,0,1)
  x3<-rnorm(nC,0,1)
  xC<-x1^2+x2^2+x3^2
  uC<-rnorm(nC,0,sd=sqrt((xC)^2))
  
  yC<-beta1+beta2*xC+uC
  olsC<-lm(yC~xC)
  sampleC<-olsC
  coefC1_100[i]<-olsC$coefficients[1]
  coefC2_100[i]<-olsC$coefficients[2]
  
  homoC<-vcov(sampleC)
  homoC
  
  varbetaC1H_100[i]<-vcov(sampleC)[1,1]
  varbetaC2H_100[i]<-vcov(sampleC)[2,2]
  varbetaC1R_100[i]<-vcovHC(sampleC,type="HC0")[1,1]
  varbetaC2R_100[i]<-vcovHC(sampleC,type="HC0")[2,2]
  
  sqC1_100<-sqrt(vcov(sampleC)[1,1])
  sqC2_100<-sqrt(vcov(sampleC)[2,2])
  vcmrC_100<- vcovHC(sampleC, type="HC0")
  vcmrC_100
  
  lowerC1_100[i] <- (coefC1_100[i]) - 1.96 * sqrt(varbetaC1H_100[i])
  upperC1_100[i] <- (coefC1_100[i]) + 1.96 * sqrt(varbetaC1H_100[i]) 
  
  lowerC2_100[i] <- (coefC2_100[i]) - 1.96 * sqrt(varbetaC2H_100[i])
  upperC2_100[i] <- (coefC2_100[i]) + 1.96 * sqrt(varbetaC2H_100[i]) 
  
  lowerC1R_100[i] <- (coefC1_100[i]) - 1.96 * sqrt(varbetaC1R_100[i])
  upperC1R_100[i] <- (coefC1_100[i]) + 1.96 * sqrt(varbetaC1R_100[i]) 
  
  lowerC2R_100[i] <- (coefC2_100[i]) - 1.96 * sqrt(varbetaC2R_100[i])
  upperC2R_100[i] <- (coefC2_100[i]) + 1.96 * sqrt(varbetaC2R_100[i]) 
  
  lowerC1_100A[i] <- (coefC1_100[i]) - qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC1H_100[i])
  upperC1_100A[i] <- (coefC1_100[i]) + qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC1H_100[i]) 
  
  lowerC2_100A[i] <- (coefC2_100[i]) - qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC2H_100[i])
  upperC2_100A[i] <- (coefC2_100[i]) + qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC2H_100[i]) 
  
  lowerC1R_100A[i] <- (coefC1_100[i]) - qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC1R_100[i])
  upperC1R_100A[i] <- (coefC1_100[i]) + qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC1R_100[i]) 
  
  lowerC2R_100A[i] <- (coefC2_100[i]) - qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC2R_100[i])
  upperC2R_100A[i] <- (coefC2_100[i]) + qt(0.975, df = sampleC$df[1]) * sqrt(varbetaC2R_100[i]) 

}

### 0.95 coverage test for z and t under homoskedasticity and robust to heteroskedasticity ###

### 0.95 covarage for beta=1 homoskedasticity NORMAL ###

CIC1_100 <- cbind(lowerC1_100, upperC1_100)
mean(CIC1_100[, 1] <= 1 & 1 <= CIC1_100[, 2])
CIC2_100 <- cbind(lowerC2_100, upperC2_100)
mean(CIC2_100[, 1] <= 1 & 1 <= CIC2_100[, 2])

### 0.95 covarage for beta=0 homoskedasticity NORMAL ###

CIC1_100zero <- cbind(lowerC1_100, upperC1_100)
mean(CIC1_100zero[, 1] <= 0 & 0 <= CIC1_100zero[, 2])
CIC2_100zero <- cbind(lowerC2_100, upperC2_100)
mean(CIC2_100zero[, 1] <= 0 & 0 <= CIC2_100zero[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty NORMAL ###

CIC1R_100 <- cbind(lowerC1R_100, upperC1R_100)
mean(CIC1R_100[, 1] <= 1 & 1 <= CIC1R_100[, 2])
CIC2R_100 <- cbind(lowerC2R_100, upperC2R_100)
mean(CIC2R_100[, 1] <= 1 & 1 <= CIC2R_100[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty NORMAL ###

CIC1R_100zero <- cbind(lowerC1R_100, upperC1R_100)
mean(CIC1R_100zero[, 1] <= 0 & 0 <= CIC1R_100zero[, 2])
CIC2R_100zero <- cbind(lowerC2R_100, upperC2R_100)
mean(CIC2R_100zero[, 1] <= 0 & 0 <= CIC2R_100zero[, 2])

### 0.95 covarage for beta=1 homoskedasticity T-STUDENT ###

CIC1_100A <- cbind(lowerC1_100A, upperC1_100A)
mean(CIC1_100A[, 1] <= 1 & 1 <= CIC1_100A[, 2])
CIC2_100A <- cbind(lowerC2_100A, upperC2_100A)
mean(CIC2_100A[, 1] <= 1 & 1 <= CIC2_100A[, 2])

### 0.95 covarage for beta=0 homoskedasticity T-STUDENT ###

CIC1_100zeroA <- cbind(lowerC1_100A, upperC1_100A)
mean(CIC1_100zeroA[, 1] <= 0 & 0 <= CIC1_100zeroA[, 2])
CIC2_100zeroA <- cbind(lowerC2_100A, upperC2_100A)
mean(CIC2_100zeroA[, 1] <= 0 & 0 <= CIC2_100zeroA[, 2])

### 0.95 covarage for beta=1 robust heteroskedasticty T-STUDENT ###

CIC1R_100A <- cbind(lowerC1R_100A, upperC1R_100A)
mean(CIC1R_100A[, 1] <= 1 & 1 <= CIC1R_100A[, 2])
CIC2R_100A <- cbind(lowerC2R_100A, upperC2R_100A)
mean(CIC2R_100A[, 1] <= 1 & 1 <= CIC2R_100A[, 2])

### 0.95 covarage for beta=0 robust heteroskedasticty T-STUDENT ###

CIC1R_100zeroA <- cbind(lowerC1R_100A, upperC1R_100A)
mean(CIC1R_100zeroA[, 1] <= 0 & 0 <= CIC1R_100zeroA[, 2])
CIC2R_100zeroA <- cbind(lowerC2R_100A, upperC2R_100A)
mean(CIC2R_100zeroA[, 1] <= 0 & 0 <= CIC2R_100zeroA[, 2])

### plot Y ###

plot(xC, yC,
     ylab="yC", xlab="xC")
abline(sampleC, col="blue")

cov(xC,uC)

### comapre variancee-covariance matrix under hoskedastic and robust to heteroskedastic ###

varC1H_100<-mean(varbetaC1H_100)
varC2H_100<-mean(varbetaC2H_100)
varC1R_100<-mean(varbetaC1R_100)
varC2R_100<-mean(varbetaC2R_100)
c(varC1H_100,varC1R_100)
c(varC2H_100,varC2R_100)

homoC
vcmrC_100

sqC1_100
sqC2_100

### t-test under beta=1 and beta=0 ###

summary(olsC)
robustolsC<-robust.summary(olsC, vcov = vcovHC(olsC), type = "HC3")

(tstatC <- (robustolsC$coefficients[2,1]-1)/robustolsC$coefficients[2,2])
(pvalueC <- (pt(tstatB2, df = 500-2, lower.tail = TRUE)*2))
(pvalueC <- (pnorm(tstatB2, lower.tail = FALSE)*2))

(tstatC0 <- (robustolsC$coefficients[2,1])/robustolsC$coefficients[2,2])
(pvalueC0 <- (pt(tstatC0, df = 500-2, lower.tail = TRUE)*2))
(pvalueC0 <- (pnorm(tstatC0, lower.tail = FALSE)*2))


#################################################################################################

################################
### FOR COMMENTS SEE PDF FILE###
################################

### compare the distribution ###

### distribution coefficients A ###

par(mfrow=c(2,2)) 

hist(coefA1_50)
hist(coefA1_500)

hist(coefA2_50)
hist(coefA2_500)

### distribution coefficients B ###

par(mfrow=c(2,2)) 

hist(coefB1_50)
hist(coefB1_500)

hist(coefB2_50)
hist(coefB2_500)

### distribution coefficients C ###

hist(coefC1_100)
hist(coefC2_100)


