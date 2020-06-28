# Edit from Fig2_Knot5.R
rm(list=ls())  
source('R_functions.R')
Rcpp::sourceCpp('Rcpp_functions.cpp')
#### simulation setting ####
p = 10
knot = 10
tol = 10^(-6)
df = 3
rate = 0.001
track = 5
M_stop = 5000
rate_GD = 0.01

dat1 = sim_data(nF = 20, p =10, seed = 1)
delta = dat1$delta
facility = dat1$facility
z = dat1$z
time = dat1$time

beta_true=matrix( rep(0, length(time)*p), nrow=length(time))
beta_true[,1]=1
beta_true[,2]=sin(3*time*pi/4)
beta_true[,3]=-1
beta_true[,4]=-(time/3)**2*exp(time/2)
beta_true[,5]=1

## Newton-Raphson method
fit_NR = TimeVarying_NR(time, delta, z, facility, knot, M_stop, rate = 0.001, tol = 10^(-6))
## GDBoost method
fit_GD = TimeVarying_GDboost(time, delta, z, facility, knot, M_stop, rate = 0.01, track = 5, tol = 10^(-6))


j = 1
ylim=range(fit_NR$beta[,j],fit_GD$beta[,j],beta_true[,j])
plot(time,fit_NR$beta[,j],type = "l",col=1,lwd=2,ylab=expression(beta),ylim = ylim)
lines(time,fit_GD$beta[,j],col=2,lwd=2)
lines(time,beta_true[,j],col=3,lwd=2)
legend("topleft",legend = c("Truth","GDBoost","NR"),lwd=c(2,2,2),col=c(3,2,1))

j = 4
ylim=range(fit_NR$beta[,j],fit_GD$beta[,j],beta_true[,j])
plot(time,fit_NR$beta[,j],type = "l",col=1,lwd=2,ylab=expression(beta),ylim = ylim)
lines(time,fit_GD$beta[,j],col=2,lwd=2)
lines(time,beta_true[,j],col=3,lwd=2)
legend("bottomleft",legend = c("Truth","GDBoost","NR"),lwd=c(2,2,2),col=c(3,2,1))


