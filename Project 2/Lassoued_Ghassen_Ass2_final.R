rm(list=ls())

#Question2.1######################################################################################################
##1)##############################################################################################################
polyroot(c(-0.3,0.4,1))
##2)see the report
##3)
###The phi_i parameters in the 'arima' function in R are positive on the right hand side of the equal sign
###so we implement 0.8 and not (-0.8) in AR part
N<-200
set.seed(123)
sim = replicate(10, arima.sim(list(order = c(1,0,2), ar = c(0.8), ma = c(0.4, -0.3)),n = 200, sd = 0.4))
matplot(sim, type = "l", ylab="X_t", xlab="t")
##4)
simulation_acf = sapply(split(sim,col(sim)), function(ts) acf(ts,plot=FALSE)$acf)
x_seq = seq(from = 0, to=length(simulation_acf[,1])-1, by = 1)
matplot(x_seq, simulation_acf, type="l",ylab="ACF", xlab="Lag", main = "ACF for realizations")
axis(1, at=x_seq)
abline( 2/sqrt(N),0, lwd=1, lty=2)
abline(-2/sqrt(N),0, lwd=1, lty=2)

##5)
simulation_pacf = sapply(split(sim,col(sim)), function(ts) pacf(ts,plot=FALSE)$acf)
x_seq_2 = seq(from = 0, to=length(simulation_pacf[,1]), by = 1)
matplot(simulation_pacf, type="l",ylab="PACF", xlab="Lag",  main = "PACF for realizations")
axis(1,at=x_seq_2)
abline(2/sqrt(N),0, lwd=1, lty=2)
abline(-2/sqrt(N),0, lwd=1, lty=2)
##6)
#sim_var = sapply(split(sim,col(sim)), function(ts) var(ts))
sim_acof = sapply(split(sim,col(sim)), function(ts) acf(ts,plot=F, type="covariance")$acf)
plot(sim_acof[1,],ylab="Variance[ X(·,w) ]", xlab = " Realization w ")
##7)
### Theoritical ACF autocorrelation function ###
rho_0<- 1
rho_1<- 0.81096
rho_2<- 0.56657
rho_list <- 3:23
rho<- c(rep(0,length(x_seq))  )
for (i in rho_list)
{
  rho[i+1]<- (0.8)^( i -2) * rho_2   # rho[i] here = rho[i-1] in manuscript
}
rho[1]<-rho_0
rho[2]<-rho_1
rho[3]<-rho_2

matplot(x_seq, simulation_acf, type="l",ylab="ACF", xlab="Lag", main = "ACF for realizations")
axis(1, at=x_seq)
abline( 2/sqrt(N),0, lwd=1, lty=2)
abline(-2/sqrt(N),0, lwd=1, lty=2)
abline(0,0, lwd=1, lty=1)
lines(x_seq,rho, type="l",lty=1, col="black")
points(x_seq,rho, cex=0.8 , col="black")
legend("topright", legend=c("Analytical ACF"),  col=c("black"), pch=c(1), pt.bg =c("black") ,lty=1, cex=0.8 )

#Question2.2##################################################################################################
A2 = read.table("/home/ghassen97/Desktop/S8/time series analysis/assignment/assignment 2/A2_sales.txt", header = TRUE)
plot(A2$Sales, ylab="Sales", xlab='time', main="A2_sales", type = "l")
points(A2$Sales, cex = .5, col = "dark blue")

#predicting 2 steps ahead:#
mu=2092

Z= A2$Sales - mu  # variable change
t=21
Z[21] <- 0.99*Z[t-1]-0.22*Z[t-2]+0.62*Z[t-4]-0.6138*Z[t-5]+0.1364*Z[t-6]
t=22
Z[22] <- 0.99*Z[t-1]-0.22*Z[t-2]+0.62*Z[t-4]-0.6138*Z[t-5]+0.1364*Z[t-6]

Y= Z+ mu
plot(Y[1:20], type = 'l', xlim = c(1,25), ylab ='time', main = 'A2_sales and 2 steps predictions',col="blue")
points(A2$Sales, col = "dark blue")
points(21, Y[21], col = 'red')
points(22, Y[22], col = 'red')
legend("topright", legend=c("observations", "predictions"), col=c("blue", "red"),pch=c(1,1), pt.bg =c("blue", "red") ,lty=1, cex=0.7)

#Confidence intervals#
k<- 1
alpha<-0.05
sigma_e<-sqrt(39508)
half_interval_t1 = qt(1-alpha/2, df=Inf ) * sigma_e*sqrt(1)
half_interval_t2 = qt(1-alpha/2, df=Inf ) * sigma_e*sqrt(1+0.99^2)
plot(Y[1:20], type = 'l', xlim = c(1,25),xlab='time', ylab ='Sales', ylim=c(1500,3000) , main = 'A2_sales and 2 steps predictions',col="blue")
points(A2$Sales, col = "dark blue")
points(21, Y[21], col = 'red')
points(21, Y[21]+half_interval_t1 , col="green" )
points(21, Y[21]-half_interval_t1 , col="green" )
points(22, Y[22], col = 'red')
points(22, Y[22]+half_interval_t2 , col="green" )
points(22, Y[22]-half_interval_t2 , col="green" )
legend("bottom", legend=c("observations", "predictions","confidence interval"), col=c("blue", "red","green"),pch=c(1,1,1), pt.bg =c("blue", "red","green") ,lty=1, cex=0.7)

#Question 2.3##################################################################################################
##1)
phi2_1<- 0.52
phi2_2<- 0.98
roots_1 <- polyroot(c(phi2_1, 1.5 , 1))
norm_roots_1<- Mod(roots_1)
roots_2 <- polyroot(c(phi2_2, 1.5 , 1))
norm_roots_2<- Mod(roots_2)
##2)
###process1###
set.seed(123)
phi2_ind<-2
break1 <- seq(0.35, 0.65, length.out = 30) # to have same number of bids, and adjust x_axis

sim1 <- replicate(100, arima.sim(list(order = c(2,0,0), ar = c(-1.5, -0.52)), n = 300, sd = 0.1)) 
sim1_arima_process1 <- sapply(split(sim1,col(sim1)), function(ts) arima(ts, order = c(2,0,0))$coef)*(-1)
sim1_arima <- sim1_arima_process1[phi2_ind,]
sim1_arima_histogram<-hist(sim1_arima,breaks=break1,plot=FALSE)
plot(sim1_arima_histogram, xlab='estimated phi2', ylab='density', ylim = c(0,20), main="Process 1" )

p1<-quantile(sim1_arima,0.975)
p2<-quantile(sim1_arima,0.025)

abline(v=p1, col='red')
abline(v=p2,col='red')

curve(dnorm(x, mean=0.52061432, sd=0.04795787),add = TRUE,col=3,lwd=3)


###process2###
break2<-break1

sim2 <- replicate(100, arima.sim(list(order = c(2,0,0), ar = c(-1.5, -0.52)), n = 300, sd = 5)) 
sim2_arima_process2 <- sapply(split(sim2,col(sim2)), function(ts) arima(ts, order = c(2,0,0))$coef)*(-1)
sim2_arima <- sim2_arima_process2[phi2_ind,]

sim2_arima_histogram<-hist(sim2_arima,breaks=break2,plot=FALSE)
plot(sim2_arima_histogram, xlab='estimated phi2', ylab='density', ylim = c(0,20),main="Process 2")
p2_1<-quantile(sim2_arima,0.975)
p2_2<-quantile(sim2_arima,0.025)
abline(v=p2_1, col='red')
abline(v=p2_2,col='red')

curve(dnorm(x, mean=0.522379960, sd=0.044206581),add = TRUE,col=3,lwd=3)


###process 3#####
break3 = seq(0.9, 1.02, length.out = 30)

sim3 <- replicate(100, arima.sim(list(order = c(2,0,0), ar = c(-1.5, -0.98)), n = 300, sd = 0.1)) 
sim3_arima_process3 <- sapply(split(sim3,col(sim3)), function(ts) arima(ts, order = c(2,0,0))$coef)*(-1)
sim3_arima <- sim3_arima_process3[phi2_ind,] # phi2

sim3_arima_histogram<-hist(sim3_arima,breaks = break3,plot=FALSE)
plot(sim3_arima_histogram, xlab='estimated phi2', ylab='density', ylim = c(0,25),main="Process 3")
p3_1<-quantile(sim3_arima,0.975)
p3_2<-quantile(sim3_arima,0.025)
abline(v=p3_1, col='red')
abline(v=p3_2,col='red')

curve(dnorm(x, mean=0.974997395, sd=0.013742083),add = TRUE,col=3,lwd=3)


###process 4#####
break4<-break3
sim4 <- replicate(100, arima.sim(list(order = c(2,0,0), ar = c(-1.5, -0.98)), n = 300, sd = 5)) 
sim4_arima_process4 <- sapply(split(sim4,col(sim4)), function(ts) arima(ts, order = c(2,0,0))$coef)*(-1)
sim4_arima <- sim4_arima_process4[phi2_ind,] # phi2

sim4_arima_histogram<-hist(sim4_arima,breaks = break4,plot=FALSE)
plot(sim4_arima_histogram, xlab='estimated phi2', ylab='density', ylim = c(0,25),main="Process 4")
p4_1<-quantile(sim4_arima,0.975)
p4_2<-quantile(sim4_arima,0.025)
abline(v=p4_1, col='red')
abline(v=p4_2,col='red')

curve(dnorm(x, mean=0.976334331, sd=0.011170152),add = TRUE,col=3,lwd=3)


#3) #4) see report 
#5)
### Process 1 ###
phi1_1_hat<-sim1_arima_process1[1,]
phi2_1_hat<-sim1_arima_process1[2,]
plot(phi1_1_hat,phi2_1_hat,ylim=c(0.35,0.65), xlim=c(1.35, 1.63), xlab="estimated phi1", ylab = "estimated phi2", main='Process 1')
points(1.5, 0.52, col="blue") #true value
legend("topleft", legend=c("true values"),  col=c("blue"), pch=c(1), pt.bg =c("blue"), cex=0.8 )

cov1<-cov(phi1_1_hat,phi2_1_hat,use = "everything",method = "pearson")
cor1<-cor(phi1_1_hat,phi2_1_hat,use = "everything",method = "pearson")

var(phi1_1_hat)
var(phi2_1_hat)

### checking stationarity for all phi1 and phi2 ###
roots_process1 <- polyroot(c(phi2_1_hat, phi1_1_hat , 1))
norm_roots_process1<- Mod(roots_process1)


### Process 2 ###
phi1_2_hat<-sim2_arima_process2[1,]
phi2_2_hat<-sim2_arima_process2[2,]
plot(phi1_2_hat,phi2_2_hat,ylim=c(0.35,0.65), xlim=c(1.35, 1.63),xlab="estimated phi1", ylab = "estimated phi2", main='Process 2')
points(1.5, 0.52, col="blue")
legend("topleft", legend=c("true values"),  col=c("blue"), pch=c(1), pt.bg =c("blue"), cex=0.8 )

cov2<-cov(phi1_2_hat,phi2_2_hat,use = "everything",method = "pearson")
cor2<-cor(phi1_2_hat,phi2_2_hat,use = "everything",method = "pearson")

var(phi1_2_hat)
var(phi2_2_hat)
### checking stationarity for all phi1 and phi2 ###
roots_process2 <- polyroot(c(phi2_2_hat, phi1_2_hat , 1))
norm_roots_process2<- Mod(roots_process2)

### Process 3 ###
phi1_3_hat<-sim3_arima_process3[1,]
phi2_3_hat<-sim3_arima_process3[2,]
plot(phi1_3_hat,phi2_3_hat,ylim=c(0.9,1.02), xlim=c(1.44, 1.54),xlab="estimated phi1", ylab = "estimated phi2", main='Process 3')
points(1.5, 0.98, col="blue")
legend("topleft", legend=c("true values"),  col=c("blue"), pch=c(1), pt.bg =c("blue"), cex=0.8 )

cov3<-cov(phi1_3_hat,phi2_3_hat,use = "everything",method = "pearson")
cor3<-cor(phi1_3_hat,phi2_3_hat,use = "everything",method = "pearson")

var(phi1_3_hat)
var(phi2_3_hat)
### checking stationarity for all phi1 and phi2 ###
roots_process3 <- polyroot(c(phi2_3_hat, phi1_3_hat , 1))
norm_roots_process3<- Mod(roots_process3)

### Process 4 ###
phi1_4_hat<-sim4_arima_process4[1,]
phi2_4_hat<-sim4_arima_process4[2,]
plot(phi1_4_hat,phi2_4_hat,ylim=c(0.9,1.02), xlim=c(1.44, 1.54),xlab="estimated phi1", ylab = "estimated phi2", main='Process 4')
points(1.5, 0.98, col="blue")
legend("topleft", legend=c("true values"),  col=c("blue"), pch=c(1), pt.bg =c("blue"), cex=0.8 )

cov4<-cov(phi1_4_hat,phi2_4_hat,use = "everything",method = "pearson")
cor4<-cor(phi1_4_hat,phi2_4_hat,use = "everything",method = "pearson")

var(phi1_4_hat)
var(phi2_4_hat)
### checking stationarity for all phi1 and phi2 ###
roots_process4 <- polyroot(c(phi2_4_hat, phi1_4_hat , 1))
norm_roots_process4<- Mod(roots_process4)


### fitting histograms with normal distribution ### to be executed before 
### the 4 lines like:  curve(dnorm(x, mean=0.976334331, sd=0.011170152),add = TRUE,col=3,lwd=3)

library("MASS")
a4<-fitdistr(phi2_4_hat, "normal")
a3<-fitdistr(phi2_3_hat, "normal")
a2<-fitdistr(phi2_2_hat, "normal")
a1<-fitdistr(phi2_1_hat, "normal")





