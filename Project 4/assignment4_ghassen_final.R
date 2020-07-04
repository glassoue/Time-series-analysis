rm(list=ls())

library(data.table)
library("dlm")
library("FKF")
library(scales)
library("plot3D")
library(plotly)
library(numDeriv)
library(MASS)
library(mvtnorm)
require(graphics)

#Question1############################################################################################################
data = read.csv(file = "/home/ghassen97/Desktop/S8/time series analysis/assignment/assignment 4/A4_Kulhuse.csv", fill = TRUE, header = TRUE)

S = data$Sal
DO = data$ODO
DT = data$DateTime
i = 1:length(S) # Time

plot(S, type = 'l', main = expression(S[t]), ylab='salinity')
points(c(1296, 4609),c(0.43,1.32), col = 'blue', pch = 10)
plot(DO, type = 'l', main = expression(DO[t]), col = 'black',ylab='dissolved oxygen')

points(c(1422), c(4.53), col = 'blue', type = 'p')
points(c(4609, 1296),c(12.60,     11.26), col = 'blue', pch = 10)

#Question2 see report ############################################################################################################

#Question3############################################################################################################

Kalman_filter = function(obs_condition, outlier_condition) {
  S = data$Sal # S_t
  N = length(S)
  X = vector(length = N) # State of salinity (state vector/number)
  X_pred = vector(length = N) # One-step predictions of X
  S_xx = vector(length = N) # This is Conditional covariance of X given all previous S
  S_yy = vector(length = N) # This is Variance of S given all previous S
  S_xx_pred = vector(length = N) #  one step prediction variance
  S_yy_pred = vector(length = N) #  one step prediction variance
  K    = vector(length = N) # Kalman gain
  outlier_ind    = vector(length = N) # Outlier index
  S_t  = vector(length = N)
  error1  = vector(length = N)
  i = 1:N # Time
  V_e = 0.01 #system variance
  V_m = 0.005 #observation variance
  
  # Initial Condition
  X[1] = S[1]
  X_pred[1] = S[1]
  S_xx[1] = V_e
  S_yy[1] = V_e + V_m
  S_xx_pred[1] = V_e
  S_yy_pred[1] = V_e + V_m
  
  # Kalman Filter
  for (t in i) {
    
    error1[t] = S[t] - X_pred[t]
    
    # fix missing observations (if obs_condition = 1) fix outliers (if outlier_condition = 1). This if statement skips the reconstruction step
    # if obs_condition == TRUE then fix missing obersvations, if outlier_condition == TRUE then skip outliers. 
    # There is no reconstruction in this if below
    if((S[t] %in% NA  && obs_condition == TRUE ) | (abs(error1[t]) > 6*sqrt(S_xx[t]) && outlier_condition == TRUE)){
      outlier_ind[t] = t
      X[t+1] = X[t]
      X_pred[t+1] = X[t+1]
      S_xx[t+1] = S_xx[t] + V_e
      S_yy[t+1] = S_xx[t+1] + V_m 
      S_xx_pred[t+1] = S_xx[t] 
      S_yy_pred[t+1] = S_xx[t+1] 
    }else{
      
      
      K[t] = S_xx[t]/S_yy[t] #Kalman gain
      
      #Reconstruction , updating
      X[t] = X_pred[t] + K[t]*(S[t] - X_pred[t]) # X_t|t
      S_xx[t] = S_xx[t] - K[t]*S_xx[t] 
      
      #Prediction
      X[t+1] = X[t]
      X_pred[t+1] = X[t+1]
      S_xx[t+1] = S_xx[t] + V_e
      S_yy[t+1] = S_xx[t+1] + V_m 
      S_xx_pred[t+1] = S_xx[t] 
      S_yy_pred[t+1] = S_xx[t+1] 
    }
  }
  
  CI_upper = X_pred + 1.96*sqrt((S_yy)) # 95% confidence interval upper 
  CI_lower = X_pred - 1.96*sqrt((S_yy)) # 95% Confidence interval lower 
  
  outlier_ind = outlier_ind[!is.na(S)]
  outlier_ind = outlier_ind[outlier_ind != 0]
  # with fkf implementation, check if Kalman_filter is ok
  mod = fkf(a0 = c(S[1]), P0 = matrix(V_e), dt = matrix(0), ct = matrix(0), Tt = matrix(1), Zt = matrix(1), 
            HHt = matrix(V_e), GGt = matrix(V_m), yt = matrix(t(S), nrow = 1, ncol = length(S)), check.input = TRUE)
  check_fkf_1 = sum(X_pred - mod$at) #check_fkf_1=5,329071e-15
  #print(check_fkf_1)
  check_fkf_2 = sum(X[i] - mod$att)#check_fkf_2=5,329071e-15
  print(check_fkf_2)
  error2 = S - X_pred[i]
  
  # Return Values
  newlist = list("X" = X, "X_pred" = X_pred, "S_xx" = S_xx, "S_yy" = S_yy, 
                 "CI_upper" = CI_upper, "CI_lower" = CI_lower, "fkf" = mod, "check_fkf_1" = check_fkf_1, 
                 "check_fkf_2" = check_fkf_2, "S_xx_pred" = S_xx_pred, "S_yy_pred" = S_yy_pred,
                 "error1" = error1, "error2" = error2, 'outlier_ind' = outlier_ind, "K" = K)
  return(newlist)
  
}
#1)######################################################
K_filter_1 = Kalman_filter(TRUE, FALSE)
plot(S, type = 'l', col = 'black', main = "One-step predictions", ylab = "Salinity")
lines(K_filter_1$X_pred, type = 'l', col = 'red')
lines(K_filter_1$CI_upper, type = 'l', col = 'blue')
lines(K_filter_1$CI_lower, type = 'l', col = 'blue')
legend(1500, 10, legend = c(expression(S[t]), expression(X[t]), "Prediction Interval"), col = c("black", "red", "blue"), lty=1:1, cex = 0.9)
#2)#################

i = 1:length(S)

S_pred_error = (S - K_filter_1$X_pred[i])/sqrt(K_filter_1$S_yy[i])
plot(S_pred_error, type = 'p', ylab = 'predictions errors', main = "standardized 1 step predictions errors",pch=16,cex=0.3)

#3)##########################

K_filter_1 = Kalman_filter(TRUE, FALSE)

plot(S, type = 'l', col = 'black', main = "One-step predictions", ylab = "Salinity",xlim = c(800,950), ylim = c(15,20))
lines(K_filter_1$X_pred, type = 'l', col = 'red')
lines(K_filter_1$CI_upper, type = 'l', col = 'blue')
lines(K_filter_1$CI_lower, type = 'l', col = 'blue')
legend("bottomright", legend = c(expression(S[t]), expression(X[t]), "Prediction Interval"), col = c("black", "red", "blue"), lty=1:1, cex = 0.9)


S_pred_error = (S - K_filter_1$X_pred[i])/sqrt(K_filter_1$S_yy[i])
plot(S_pred_error, type = 'p', col = 'black', xlim = c(800,950), ylim = c(-20,20), ylab = 'predictions errors', main = "standardized 1 step predictions errors",pch=16,cex=0.4)

#4)#########################
#using function to get X
K_filter_1$X[5000] 
K_filter_1$S_xx[5000]

# using library FKF
K_filter_1$fkf$att[5000]  #same values
K_filter_1$fkf$Ptt[5000]

#Question4############################################################################################################
#1)#################
K_filter_no_outlier = Kalman_filter(TRUE, TRUE)

plot(S, type = 'l', col = 'black', main = "One-step predictions", ylab = "Salinity",xlim = c(800,950), ylim = c(15,20))
lines(K_filter_no_outlier$X_pred, type = 'l', col = 'red')
lines(K_filter_no_outlier$CI_upper, type = 'l', col = 'blue')
lines(K_filter_no_outlier$CI_lower, type = 'l', col = 'blue')
legend("bottomright", legend = c(expression(S[t]), expression(X[t]), "Prediction Interval"), col = c("black", "red", "blue"), lty=1:1, cex = 0.9)

S_pred_error = (S - K_filter_no_outlier$X_pred[i])/sqrt(K_filter_no_outlier$S_yy[i])
plot(S_pred_error, type = 'p', col = 'black', xlim = c(800,950), ylim = c(-20,20), ylab = 'predictions errors', main = "standardized 1 step predictions errors",pch=16,cex=0.4)

#2)################
K_filter_no_outlier$outlier_ind[1:5] #Outliers from the output of tje kalman filter K_filter_no_outlier = Kalman_filter(1,0,1)
#3)################
length(K_filter_no_outlier$outlier_ind) #Number of skipped observations excluding the missing observations (NA)
#4)###############
K_filter_no_outlier$X[5000] 
K_filter_no_outlier$S_xx[5000] 

#Question5############################################################################################################

#Question6############################################################################################################
#see report





