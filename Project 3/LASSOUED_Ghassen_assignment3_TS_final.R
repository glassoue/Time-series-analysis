rm(list=ls())
library(forecast)
library(car)

#Question1####################################################################################################
data = read.table("/home/ghassen97/Desktop/S8/time series analysis/assignment/assignment 3/2020_A3_co2.txt", header = TRUE)
co2 = data$co2
t   = data$time
#test with last 24 observations, length(t)=601

plot(t[1:576], co2[1:576],xlim = c(1970,2020), type='l', main = expression(paste(Y[t], "  and  ", Y_test[t])), xlab = expression(t), col = 'black')
lines(t[577:601], co2[577:601], type='l', col='red')
Y = co2[1:576]
legend("topleft", legend = c(expression(Y[t]), expression(Y_test[t])), col = c('black', 'red'),  lty=1:1)

#Question2###################################################################################################
#attach(mtcars)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(Y, type = 'l', main = expression(Y[t]))
acf(Y,lag.max = 50, main = expression(SACF(Y[t])))
pacf(Y,lag.max = 50,  main = expression(SPACF(Y[t])))

#differencing Y_t
p= 1
H = diff(Y, lag = p, differences = 1)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(H, type = 'l')
acf(H,lag.max = 50, main = expression(SACF(H[t])))
pacf(H,lag.max = 50,  main = expression(SPACF(H[t])))

#differencing H_t withperiod p=12
p = 12
Z = diff(H, lag = p, differences = 1)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(Z, type = 'l')
acf = acf(Z,lag.max=50,main = expression(SACF(Z[t])))
pacf =pacf(Z,lag.max=50, main = expression(SPACF(Z[t])))

p = 1
X = diff(Z, lag = p, differences = 1)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(X, type = 'l')
acf = acf(X, lag.max = 50, main = expression(SACF(X[t])))
pacf =pacf(X, lag.max = 50, main = expression(SPACF(X[t])))

#Question3###################################################################################################
#see report
#Question4###################################################################################################
##definition of function that will be used######

#from R example week7: 
my.tsdiag = function(model, nlag){  
  oldpar = par(mfrow=c(4,1), mgp=c(2,0.7,0), mar=c(3,3,1.5,1))
  on.exit(par(oldpar))
  plot(fitted(model), col = 'blue', type = 'l')
  lines(model$x, col = 'red', type = 'l')
  model = model$residuals
  acf(model, lag.max = 50)
  pacf(model, lag.max = 50)
  
  pval = sapply(1:nlag, function(i) Box.test(model, i, type = "Ljung-Box")$p.value)
  plot(1L:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1), main = "p values for Ljung-Box statistic")
  abline(h = 0.05, lty = 2, col = "blue")
  
  return(pval)
}

# Test whether parameter=0
test_parameter_zero = function(model, a){
  TS = model$coef/(sqrt(abs(diag(model$var.coef))))
  pval = 2*(1-pt(abs(TS), df = (length(model$residuals) - length(model$coef)-1)))
  
  decision = symnum(pval[1:(length(pval))], decision(0, a, Inf), decision("Reject", "Accept"))
  print(pval[1:(length(pval))])
  print(decision)
}

# Stopping criteria function
stopping_criteria_final = function(model, nlag){
  pval = my.tsdiag(model, nlag)
  pval_lagk = sum((pval < 0.05))
  if(pval_lagk > 0)
  {
    I = 500 * pval_lagk
  }else{
    I = 0
  }
  ni = length(model$coef) - 1
  #stopping_criteria = model$aic + ni + I , this is the remark in report in Q3.3 regardind ni
  stopping_criteria = model$aic  + I
  print(I)
  return(stopping_criteria)
}


###
t = 1:(576)
#model_1: 
model1 = Arima(Y[1:(576)], order = c(1,1,1), seasonal = list(order=c(0,1,1), period = 12), xreg = t) 
#model_2
model2 = Arima(Y[1:(576)], order = c(2,1,1), seasonal = list(order=c(0,1,1), period = 12), xreg = t)

#model_3
model3 = Arima(Y[1:(576)], order = c(1,2,1), seasonal = list(order=c(0,1,1), period = 12), xreg = t)

#model_4: 
model4 = Arima(Y[1:(576)], order = c(1,2,2), seasonal = list(order=c(0,1,1), period = 12), xreg = t) 

#model_5= Model_4 with xreg =data$time[t]
model5 = Arima(Y[1:(576)], order = c(1,2,2), seasonal = list(order=c(0,1,1), period = 12), xreg = data$time[t])

#### test whether parameter=0 using test_parameter_zero
test_parameter_zero(model1,1)
test_parameter_zero(model2,1)
test_parameter_zero(model3,1)
test_parameter_zero(model4,1)
test_parameter_zero(model5,1)
####Ljung-box test using my.tsdiag
my.tsdiag(model1,50)
my.tsdiag(model2,50)
my.tsdiag(model3,50) #fail, has several p_values < 0.05
my.tsdiag(model4,50) 
my.tsdiag(model5,50)
###stopping_criteria_final
stopping_criteria_final(model1,50)
stopping_criteria_final(model2,50)
stopping_criteria_final(model4,50)
stopping_criteria_final(model5,50)


nlag = 50 
#model1
qqPlot(model1$residuals, main = "model1", ylab ="Sample Quantiles",xlab = "Theoretical Quantiles")
dat = model1
dat = dat$residuals
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
acf(dat, lag.max = 50, main = "model1")
pacf(dat, lag.max = 50, main = "model1")
pval = sapply(1:nlag, function(i) Box.test(dat, i, type = "Ljung-Box")$p.value) #from week7 R example
plot(1L:nlag, pval, xlab = "lag", ylab = "p_value", ylim = c(0,1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")


#model2
qqPlot(model2$residuals, main = "model2", ylab ="Sample Quantiles",xlab = "Theoretical Quantiles")
dat = model2
dat = dat$residuals
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
acf(dat, lag.max = 50, main = "model2")
pacf(dat, lag.max = 50, main = "model2")
pval = sapply(1:nlag, function(i) Box.test(dat, i, type = "Ljung-Box")$p.value) #from week7 R example
plot(1L:nlag, pval, xlab = "lag", ylab = "p_value", ylim = c(0,1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")


#model4
qqPlot(model4$residuals, main = "model4",ylab ="Sample Quantiles",xlab = "Theoretical Quantiles")
dat = model4
dat = dat$residuals
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
acf(dat, lag.max = 50, main = "model4")
pacf(dat, lag.max = 50, main = "model4")
pval = sapply(1:nlag, function(i) Box.test(dat, i, type = "Ljung-Box")$p.value) #from week7 R example
plot(1L:nlag, pval, xlab = "lag", ylab = "p_value", ylim = c(0,1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")

#model5
qqPlot(model5$residuals, main = "model5", ylab ="Sample Quantiles",xlab = "Theoretical Quantiles")
dat = model5
dat = dat$residuals
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
acf(dat, lag.max = 50, main = "model5")
pacf(dat, lag.max = 50, main = "model5")
pval = sapply(1:nlag, function(i) Box.test(dat, i, type = "Ljung-Box")$p.value) #from week7 R example
plot(1L:nlag, pval, xlab = "lag", ylab = "p_value", ylim = c(0,1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")


#Question5###################################################################################################

#model1
new_t= 577:601
forecast_model1<-forecast(model1,h=25,level=0.95, xreg = new_t)
lower_model1<-forecast_model1$lower
upper_model1<-forecast_model1$upper
plot(forecast_model1, main="model 1", xlim= c(550, 610),ylim=c(360,420),xlab = 't',ylab="co2")
lines(lower_model1, type = 'l', col='red' )
lines(upper_model1, type = 'l', col='red' )
lines(1:601, data$co2,type ='l',col='black')
legend("bottomleft", legend = c("prediction on 24 test observation","confidence interval","observations"), col = c('blue', 'red','black'),  lty=1:1)

predictions_model1<-forecast_model1$mean
#SSE 
SSE_model1 <- sum ( (data$co2[577:601] - predictions_model1)^2 )

#model5
xreg_model5<-data$time[new_t]
forecast_model5<-forecast(model5,h=25,level=0.95, xreg = xreg_model5)
lower_model5<-forecast_model5$lower
upper_model5<-forecast_model5$upper
plot(forecast_model1, main="model 5",xlim= c(550, 610),ylim=c(360,420),xlab = 't',ylab="co2")
lines(lower_model5, type = 'l', col='red' )
lines(upper_model5, type = 'l', col='red' )
lines(1:601, data$co2,type ='l',col='black')
legend("bottomleft", legend = c("prediction on 24 test observation","confidence interval","observations"), col = c('blue', 'red','black'),  lty=1:1)

predictions_model5<-forecast_model5$mean
#SSE 
SSE_model5 <- sum ( (data$co2[577:601] - predictions_model5)^2 )



#Question6###################################################################################################

#model1
long<-844
new_t_long= 577:long
h_value<- long -577 + 1
forecast_model1_long<-forecast(model1,h=h_value,level=0.95, xreg = new_t_long)
# the first time model1 reaches 460 is in t=844
predictions_model1_long <- forecast_model1_long$mean
ind<-0
for (i in 1:length(predictions_model1_long) )
{
  if (predictions_model1_long[i]>460)
  {
    ind<-ind+1  #ind=1 making sure no 460 was missed
  }
} 

lower_model1_long<-forecast_model1_long$lower
upper_model1_long<-forecast_model1_long$upper
plot(forecast_model1_long, main="model 1",xlim=c(450,900), xlab = 't',ylab="co2")
lines(lower_model1_long, type = 'l', col='red' )
lines(upper_model1_long, type = 'l', col='red' )
legend("bottomright", legend = c("prediction","confidence interval","observations"), col = c('blue', 'red','black'),  lty=1:1)



#model5
# first step, define the regressor 
h_value=220
t = 577:(577+h_value-1) #t=577:796
delta_time = data$time[2:576] - data$time[1:575] # Delta
delta_time = delta_time[1:12]
delta_period_12 = (1:12)*delta_time
delta_period_12 = rep(delta_period_12, ceiling(h_value/12+1) ) 
delta_period_12 = delta_period_12[3:length(delta_period_12)]

year_start = 2018
xreg_model5_long = 0 
for (i in seq(from = year_start, to= ceiling(2018+h_value/12), by = 1)) {
  years = rep(year_start, 12)
  xreg_model5_long = append(xreg_model5_long, years, length(xreg_model5_long))
  year_start = year_start + 1
}
xreg_model5_long = xreg_model5_long[4:length(xreg_model5_long)]
xreg_model5_long = xreg_model5_long + delta_period_12
xreg_model5_long = xreg_model5_long[1:h_value]



h_value<-h_value

forecast_model5_long<-forecast(model5,h=h_value,level=0.95, xreg = xreg_model5_long)

predictions_model5_long <- forecast_model5_long$mean

ind<-0
for (i in 1:length(predictions_model5_long) )
{
  if (predictions_model5_long[i]>460)
  {
    ind<-ind+1  #indicator to make sure that 460 was reached once
  }
} 


lower_model5_long<-forecast_model5_long$lower
upper_model5_long<-forecast_model5_long$upper
plot(forecast_model5_long, main="model 5",xlim=c(450,900), xlab = 't',ylab="co2")
lines(lower_model5_long, type = 'l', col='red' )
lines(upper_model5_long, type = 'l', col='red' )
#lines(1:601, data$co2,type ='l',col='black')
legend("bottomright", legend = c("prediction","confidence interval","observations"), col = c('blue', 'red','black'),  lty=1:1)

















