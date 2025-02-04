rm(list=ls())

### ADVANCED ECONOMETRICS ###
### HOMEWORK 1 ###
### Fabio Marcaurelio (599027)###


################################
### FOR COMMENTS SEE PDF FILE###
################################


### 2.1- Plot the relationship between Y and X. Is there a linear relationship between the two variables? ### 


set.seed(599027)
n=1000
X<-rnorm(n,0,1)
Z<-rnorm(n,0,1)
Y<-(X^2 + 2*Z)
plot(X,Y,main = "relationship between Y and X", xlab = "X", ylab = "Y")
plot(X,Z,main = "relationship between X and Z", xlab = "X", ylab = "Z")

### 2.2- Evaluated E(Y) and E(XY) using the simulated draws ###

mean(Y)
mean(Y*X)

### 2.3- Evaluate cov(Y, X) using the simulated draws ###

cov(X,Y)

### 3- Use R to generate n = 100 observations from the following model (i = 1, ..., 100) ###

set.seed(599027)  
n1<- 100    
beta_1<- 1
beta_2<- 0.3
xi<-rnorm(n1,mean=3,sd=2)
u1<-rnorm(n1,0,1)
u2<-rnorm(n1,0,1)
u3<-rnorm(n1,0,1)
ui<-u1^2+u2^2+u3^2-3
yi<-beta_1+beta_2*xi+ui
yi
mean(yi)
plot(yi)
plot(xi)
mean(ui)
cov(xi,ui)

### 4- Using the data generated in the previous model to estimate β1 and β2 using OLS. Do not use the built-in command lm, but rather compute the coefficients using the formulas developed during classes.###

cov(yi,xi)
var(xi)
mean(yi)
mean(xi)
beta2hat<-cov(yi,xi)/var(xi)
beta1hat<-mean(yi)-beta2hat*mean(xi)
cbind(beta1hat,beta2hat)

### 5- Check your results using the command lm ###

ols<-lm(yi~xi)
ols      
library(summariser)
library(stargazer)
summary(ols)
stargazer(ols)
cbind(beta1hat,beta2hat)
plot(ols)

### 6 - Since you generated the data yourself, you know which assumptions hold for the model above. Answer the following questions providing a brief (1 line max) justification ###

mean(ols$residuals) 
mean(ols$residuals^2)
mean(ui^2)
var(ols$residuals)
var(ui)

### 7- ###

set.seed(599027)
V<-as.matrix(rnorm(2, mean=0,sd=1),ncol(1))
V
row1<-c(2,-1)
row2<-c(-1,3)
sigma<-rbind(row1,row2)
sigma1<-chol(sigma)
sigma1
t(sigma1)%*%sigma1
n=1000
beta1<-2
beta2<-2
x1<-rnorm(n,0,1)
x2<-rnorm(n,0,1)
u<-rnorm(n,0,1)

y<-beta1*x1+beta2*x2+u
mean(y)
ols<-lm(y~x1+x2)

vcov(ols)
chol(vcov(ols))
t(vcov(ols))%*%chol(vcov(ols))
