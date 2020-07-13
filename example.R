# Edit from Fig2_Knot5.R
rm(list=ls())  

source('R_functions.R')
Rcpp::sourceCpp('Rcpp_functions.cpp')
#### simulation setting ####
p = 5
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

## Newton-Raphson method: it may run out of memory, so this part of code is commented out
t0 = proc.time()
# fit_NR = NULL
# knot_set=quantile(time[delta],prob=seq(1:(knot-4))/(knot-3))
# library(survival)
# fit_NR = coxph(Surv(time, delta)~ tt(z[,1])+tt(z[,2])+tt(z[,3])+tt(z[,4])+tt(z[,5]),
#                        tt=function(z,time,...) {mtrx = model.matrix(~z)[,-1]
#                        t=time
#                        knots= knot_set
#                        X = matrix(splines::bs(time, df = knot, knot = knot_set, intercept = TRUE, degree = 3),nrow=length(time))
#                        mtrx * X
#                        },eps =tol)
t1 = proc.time()
t1 - t0
## GDBoost method
fit_GD = TimeVarying_GDboost(time, delta, z, facility, knot, M_stop = 500, rate = 0.05, tol = 10^(-5))
t2 = proc.time()
t2 - t1

j = 1
ylim=range(fit_GD$beta[,j],beta_true[,j])+c(-0.5,0.5)
plot(time,beta_true[,j],type = "l",col=1,lwd=2,ylab=expression(beta),ylim = ylim)
lines(time,fit_GD$beta[,j],col=2,lwd=2)
#lines(time,fit_NR$beta[,j],col=3,lwd=2)
#legend("topleft",legend = c("Truth","GDBoost","NR"),lwd=c(2,2,2),col=c(1,2,3))

j = 4
ylim=range(fit_GD$beta[,j],beta_true[,j])+c(-0.5,0.5)
plot(time,beta_true[,j],type = "l",col=1,lwd=2,ylab=expression(beta),ylim = ylim)
lines(time,fit_GD$beta[,j],col=2,lwd=2)
#lines(time,fit_NR$beta[,j],col=3,lwd=2)
#legend("bottomleft",legend = c("Truth","GDBoost","NR"),lwd=c(2,2,2),col=c(1,2,3))





#fit_NR = TimeVarying_NR(time, delta, z, facility, knot, M_stop, rate = 0.001, tol = 10^(-6))
