rm(list=ls())

## reading data and regarding the first line in the file as names:
df <- read.table("/home/ghassen97/Desktop/S8/time series analysis/assignment 1/A1_co2.txt", header=TRUE)
X<-df[,1:3]
Y<-df[,4] #Y=co2
Y=data.matrix(Y)
#Q1
## plotting
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[,3], Y, xlab='Time', ylab='co2')

## Testing with last 25 observations
y_train<-Y[1:576]
n<-length(y_train)
p<-12 #period in months
p_course=4 #number of parameters

Xols <- matrix(0,nrow=n,ncol=4)
ft <- function(i) rbind( 1, X[i,3] , sin(((2*pi)/p)* X[i,3]), cos (((2*pi)/p)* X[i,3]) )
for (i in 1:n) {
  Xols[i,]<- ft(i)
} 

#Q2 ########################################################################################
## Finding parameters manually: 
theta_hat <- solve(t(Xols)%*%Xols, t(Xols)%*%y_train) 
SIGMA=diag(n)   #MATRIX SIGMA

##Q2.2 uncertainty on theta ########################################################################################
y_hat_ols= Xols %*% theta_hat
residuals= y_train - y_hat_ols 
theta_uncertainty<- solve( t(Xols)%*%Xols, t(Xols)%*%residuals ) # theta_hat - theta

sigma2_hat= t(residuals) %*% residuals / (n-p_course)
sigma_hat=sqrt(sigma2_hat)

#plot the fitted values together with data
#todo
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[,3], Y, xlab='Time', ylab='co2', main="OLS ")
lines(X[1:n ,3], y_hat_ols, col='red' ,type ="l",lwd=2)
legend("topleft", legend=c("y", "y_hat_ols"), col=c("black", "red"), lty=1, cex=1)

## Checking distribution of residuals:
par(mfrow=c(1,1),mar=c(3,3,2,1),mgp=c(2,0.7,0))
hist(residuals, probability=T, breaks=40, main="Distribution of residuals in OLS")
curve(dnorm(x, mean=0, sd=sigma_hat),add = TRUE,col=3,lwd=3)

#Q2.4 ########################################################################################
#Q2.5 :  relaxation algorithm ########################################################################################
SIGMA=diag(n)
ind<-0
for (i in 1:5){
  
  theta_hat_relax <- solve(t(Xols) %*% solve(SIGMA) %*% Xols) %*% t(Xols)%*%solve(SIGMA)%*%y_train 
  y_hat_relax<- Xols %*% theta_hat_relax
  residuals_relax<- y_train - y_hat_relax
  rho= cor(residuals_relax[2:n], residuals_relax[1:n-1] )
  
  #SIGMA definition
  for (m in 1:n){
    for (l in 1:n){
      SIGMA[m,l]= rho ^ abs(m-l)
    }
  }
}
#Q2.6########################################################################################
theta_hat_wls <- solve(t(Xols) %*% solve(SIGMA) %*% Xols) %*% t(Xols)%*%solve(SIGMA)%*%y_train
y_hat_wls= Xols %*% theta_hat_wls
residuals_wls= y_train - y_hat_wls 
theta_wls_uncertainty<- solve( t(Xols)%*%Xols, t(Xols)%*%residuals_wls )

sigma2_hat_wls= t(residuals_wls)%*%SIGMA %*% residuals_wls / (n-p_course) #3.44 p 39
sigma_hat_wls=sqrt(sigma2_hat_wls)


#Q2.7 
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[,3], Y, xlab='Time', ylab='co2', main="WLS ")
lines(X[1:n ,3], y_hat_wls, col='red' ,type ="l",lwd=2)
legend("topleft", legend=c("y", "y_hat_wls"), col=c("black", "red"), lty=1, cex=1)

## Checking distribution of residuals:
par(mfrow=c(1,1),mar=c(3,3,2,1),mgp=c(2,0.7,0))
hist(residuals_wls, probability=T, breaks=40, main="Distribution of residuals in WLS")
curve(dnorm(x, mean=0, sd=sigma_hat_wls),add = TRUE,col=3,lwd=3)

##########################################################################################

##LOCAL LINEAR TREND
period <- 1970.125 - 1970.042 #l in course
#Q3.1
#j = t
lambda<- 0.9
L <- matrix(c(1,1,0,0,0,1,0,0,0,0,cos(2*pi/p),-sin(2*pi/p),0,0,sin(2*pi/p),cos(2*pi/p)),ncol=4)
Linv <- solve(L)
f <- function(j) rbind( 1, j, sin((2*pi/p)*j), cos((2*pi/p)*j) )
#Q3.2
#we are going to skip 10 obs after, so we'll initialize with init = 10 : #week3_exercise_sol_R

init <- 10
## FNinit & hNinit (First 10 observations)
F_n <- matrix(0, nrow=4, ncol=4)
h_n <- matrix(0,nrow=4,ncol=1)
for ( j in 0: (init-1) )
{
  F_n<- F_n + ( lambda^j *f(-j)%*%t(f(-j)) )
  h_n<- h_n + ( lambda^j *f(-j) * y_train[init-j]  )
}

## Allocating space
theta.all <- matrix(NA,ncol=p_course, nrow=n)
sigma.all <- matrix(NA,ncol=1, nrow=n)
epsilon<- matrix(NA,ncol=1, nrow=n)
sd.theta.all <- matrix(NA,ncol=p_course, nrow=n)
y_train_hat_n <- matrix(NA,ncol=1,nrow=n)
predict_int <- numeric(n)
alpha<-0.05

#T= (1-lambda^n)/(1-lambda)

## Solving at time init
theta.hat <- solve(F_n, h_n)
theta.all[init,] <- solve(F_n, h_n) #theta10

y_test=matrix(NA,ncol=1,nrow=length(Y))

## Looping over the remaining observations
for (i in (init+1):n){
  F_n <- F_n + (lambda^(i-1)) * f(- (i-1)) %*% t(f(- (i-1)))
  h_n <- (lambda) * Linv %*% h_n + f(0)*Y[i]
  theta.hat <-  solve(F_n, h_n)  
  theta.all[i,] <- theta.hat
  
  ## Adding uncertainty information
  y_train_hat_n[i]<- t(f(period)) %*% theta.hat 
  yhat=t(f(period)) %*% theta.hat
  #  y_train_hat_n[1:i]<-cbind( 1, (-i+1):0  , sin((2*pi/p)*((-i+1):0)), cos((2*pi/p)*((-i+1):0)) ) %*% theta.hat #false
  
  sigma_square_hat_upper <- 0
  T <- 0
  for(j in 1:(i-1))
  {
    sigma_square_hat_upper= sigma_square_hat_upper + lambda^(i-1-j)*(y_train[j]-yhat)^2
    T <- T + lambda^(i-1-j)
  }
  sigma_square_hat <- sigma_square_hat_upper/(T-p_course)
  epsilon[i] <- y_train[i] - y_train_hat_n[i]
  sigma.all[i] <- sqrt(sigma_square_hat)
  ## Estimating s.d. of estimated parameters
  sd.theta.all[i,] <- sigma.all[i] * sqrt(diag(solve(F_n))) 
  predict_int[i] <- qt( 1 - alpha/2, (i-1) -p_course  ) * sigma.all[i-1] * sqrt( 1+ t(f(period))%*% solve(F_n) %*% f(period) ) 
  
}

#Q3.3 
plot( X[1:n,3] ,epsilon , main="prediction errors = epsilon", xlab="time ", ylab="prediction errors ", pch=19)
plot( X[1:n,3] ,sigma.all , main="sigma ", xlab="time ", ylab="sigma ", pch=19)
#Q3.4
## 95% prediction interval 
y_train_hat_n_low <- y_train_hat_n - predict_int
y_train_hat_n_up  <- y_train_hat_n + predict_int

par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[,3], Y, xlab='Time', ylab='co2', main="95% prediction ")
lines(X[1:n ,3], y_train_hat_n, col='red' ,type ="l",lwd=2)
lines(X[1:n ,3], y_train_hat_n_up, col ='blue',lwd=2 )
lines(X[1:n ,3], y_train_hat_n_low, col ='green',lwd=2 )
legend("topleft", legend=c("y_up", "y","y_hat","y_low"), col=c("blue","black", "red","green"), lty=1, cex=1)

#Q3.5 ######################################## TEST
## Allocating space
theta_test.all <- matrix(NA,ncol=p_course, nrow=length(Y))
sigma_test.all <- matrix(NA,ncol=1, nrow=length(Y))
epsilon_test<- matrix(NA,ncol=1, nrow=length(Y))
sd.theta_test.all <- matrix(NA,ncol=p_course, nrow=length(Y))
y_test <- matrix(NA,ncol=1,nrow=length(Y))
predict_test_int <- numeric( length(Y) )
alpha<-0.05

## Looping over the test observations
for (i in (n+1):length(Y) ){
  F_n <- F_n + (lambda^(i-1)) * f(- (i-1)) %*% t(f(- (i-1)))
  h_n <- (lambda) * Linv %*% h_n + f(0)*Y[i]
  theta.hat <-  solve(F_n, h_n)  
  theta_test.all[i,] <- theta.hat
  
  ## Adding uncertainty information
  y_test[i]<- t(f(period)) %*% theta.hat 
  yhat=t(f(period)) %*% theta.hat
  sigma_square_hat_upper <- 0
  T <- 0
  for(j in 1:(i-1))
  {
    sigma_square_hat_upper= sigma_square_hat_upper + lambda^(i-1-j)*(Y[j]-yhat)^2
    T <- T + lambda^(i-1-j)
  }
  sigma_square_hat <- sigma_square_hat_upper/(T-p_course)
  epsilon_test[i] <- Y[i] - y_test[i]
  sigma_test.all[i] <- sqrt(sigma_square_hat)

  ## Estimating s.d. of estimated parameters
  sd.theta_test.all[i,] <- sigma.all[i] * sqrt(diag(solve(F_n))) 
  predict_test_int[i] <- qt( 1 - alpha/2, (i-1) -p_course  ) * sigma_test.all[i-1] * sqrt( 1+ t(f(period))%*% solve(F_n) %*% f(period) ) 
  
}
i<- (n+1)
F_n <- F_n + (lambda^(i-1)) * f(- (i-1)) %*% t(f(- (i-1)))
h_n <- (lambda) * Linv %*% h_n + f(0)*Y[i]
qt( 1 - alpha/2, (i-1) -p_course  ) * sigma_test.all[i-1] * sqrt( 1+ t(f(period))%*% solve(F_n) %*% f(period) )



## 95% prediction interval from 2010 
y_test_low <- y_test - predict_test_int
y_test_up  <- y_test + predict_test_int

###searching for index of 2010 
index<- 0
ind <- FALSE
while (ind == FALSE ) 
{
  index<-index+1
  ind <- X[index,1] == 2010
}

nn <- length(Y)
#Train since 2010
index<-420
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[index:nn,3],xlim=c(2009,2020), Y[index:nn],ylim=c(370,435), xlab='Time', ylab='co2', main="95% prediction of training data and test data ")
lines(X[index:n ,3], y_train_hat_n[index:n], col='red' ,type ="l",lwd=2)
lines(X[index:n ,3], y_train_hat_n_up[index:n], col ='green',lwd=2 )
lines(X[index:n ,3], y_train_hat_n_low[index:n], col ='green',lwd=2 )
legend("topleft", legend=c("y_up", "y","y_hat","y_low"), col=c("green","black", "red","green"), lty=1, cex=0.8)

#Test
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
lines(X[(n+1):nn ,3], y_test[(n+1):nn] , col='orange' ,type ="l",lwd=2)
lines(X[(n+1):nn ,3], y_test_up[(n+1):nn], col ='blue',lwd=2 )
lines(X[(n+1):nn ,3], y_test_low[(n+1):nn], col ='blue',lwd=2 )
legend("bottomright", legend=c("y_up_test", "y_test","y_hat_test","y_low_test"), col=c("blue","black", "yellow","blue"), lty=1, cex=0.8)

##Q3.6###################################################################"
# predicted        ; data          ; confidence interval
y_test[nn - 25 +1] ; Y[nn - 25 +1] ; predict_test_int[nn - 25 +1]   # 1 month later 
y_test[nn - 25 +2] ; Y[nn - 25 +2] ; predict_test_int[nn - 25 +2]   # 2 months later
y_test[nn - 25 +6] ; Y[nn - 25 +6] ; predict_test_int[nn - 25 +6]   # 6 months later
y_test[nn - 25 +12]; Y[nn - 25 +12]; predict_test_int[nn - 25 +12]  # 12months later
y_test[nn - 25 +24]; Y[nn - 25 +24]; predict_test_int[nn - 25 +24]  # 24 months later 

#Q3.7########################"
# compare with test data and comment 

#Q3.8######################### 
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(X[1:nn,3],xlim=c(2009,2020), Y[1:nn],ylim=c(370,450), xlab='Time', ylab='co2', main="Data & estimated mean")
lines(X[1:n ,3], theta.all[1:n], col='red' ,type ="l",lwd=2)
lines(X[(n+1):nn,3] , theta_test.all[(n+1):nn], col='yellow',type ="l",lwd=2)
legend("topleft", legend=c("y", "alpha_model_train","alpha_model_test"), col=c("black","red", "yellow"), lty=1, cex=0.8)

##########################################################################################
#Q4.1
burning_period <-100
length_alpha<-100
alpha_seq <- seq(from = 0.1, to = 0.99, length = length_alpha)
SSE_seq<- matrix(0,ncol=1,nrow=length_alpha)
#SSE_seq<- numeric(length_alpha)
for (k in 1:length_alpha )
{
  lambda <- 1 - alpha_seq[k]
  for (i in (init+1):(burning_period+init)){
    F_n <- F_n + (lambda^(i-1)) * f(- (i-1)) %*% t(f(- (i-1)))
    h_n <- (lambda) * Linv %*% h_n + f(0)*Y[i]
    theta.hat <-  solve(F_n, h_n)  
    theta.all[i,] <- theta.hat
    
    ## Adding uncertainty information
    y_train_hat_n[i]<- t(f(period)) %*% theta.hat 
    yhat=t(f(period)) %*% theta.hat
    #  y_train_hat_n[1:i]<-cbind( 1, (-i+1):0  , sin((2*pi/p)*((-i+1):0)), cos((2*pi/p)*((-i+1):0)) ) %*% theta.hat #false
    
    sigma_square_hat_upper <- 0
    T <- 0
    for(j in 1:(i-1))
    {
      sigma_square_hat_upper= sigma_square_hat_upper + lambda^(i-1-j)*(y_train[j]-yhat)^2
      T <- T + lambda^(i-1-j)
    }
    sigma_square_hat <- sigma_square_hat_upper/(T-p_course)
    epsilon[i] <- y_train[i] - y_train_hat_n[i]
    #sigma.all[i] <- sqrt(sigma_square_hat)
    ## Estimating s.d. of estimated parameters
    #sd.theta.all[i,] <- sigma.all[i] * sqrt(diag(solve(F_n))) #also for global to change?
    #predict_int[i] <- qt( 1 - alpha/2, (i-1) -p_course  ) * sigma.all[i-1] * sqrt( 1+ t(f(period))%*% solve(F_n) %*% f(period) ) 
  }
  S<-0
  for (w in (init+1):length(epsilon))
  {
    S<-S+ (epsilon[w])^2
  }
    
  SSE_seq[k] <- S
}
par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(alpha_seq,SSE_seq, xlab="alpha", ylab="SSE", main="SSE = f(alpha)" )
#lines(X[1:n ,3], theta.all[1:n], col='red' ,type ="l",lwd=2) 




