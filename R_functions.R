library(survival)
library(mvtnorm)

# Construct AR1 covariance matrix
AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}

# stratified sampling
stra_sampling<-function(size, prop){
  m<-length(size)
  mark<-cumsum(size)
  total<-sum(size)
  sample.id<-NULL
  sample.id<-c(sample.id,sample.int(size[1],size=floor(size[1]/prop),replace=FALSE))
  for(i in 2:m){
    sample.id<-c(sample.id,mark[i-1]+sample.int(size[i],size=floor(size[i]/prop),replace=FALSE))
  }
  return(sample.id)
}


# data simulation 
sim_data <- function(nF = 20, p =10, seed = 1){
  set.seed(seed)
  n_f = rpois(nF, lambda = 500)  #sample size for each facility
  N=sum(n_f)
  gamma = rnorm(nF, mean=0, sd=0.5)
  range(gamma)
  gamma_subject=rep(gamma,n_f)
  F_pre=1:nF
  facility=rep(F_pre, n_f)
  ############generate data########################
  Sigma_z1<-AR1(0.6,p)
  z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  z_012_rare=function(x){
    # U=runif(2, 0.97,1)
    U=runif(2, 0.9,0.975)
    U=U[order(U)]
    x2=quantile(x,prob=U)
    x3=x
    x3[x<x2[1]]=0
    x3[x>x2[2]]=2
    x_index=(((x<x2[1])+(x>x2[2]))<1) 
    x3[x_index]=1
    return(x3)
  }
  z[,1:3]=apply(z[,1:3],2, z_012_rare)
  MAF=apply(z,2,mean)
  range(MAF) #0.0769688 0.1778108
  U=runif(N, 0,1)
  pre_time=rep(0, N)
  for (i in 1:(N)) {
    f=function(t) {
      integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]-(x/3)**2*exp(x/2)*(x<3)*z[i,4]+z[i,5])}
      Lambda=integrate(integrand, lower = 0, upper = t)$value
      Lambda+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      pre_time[i]=4
    }
    else pre_time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  pre_censoring=runif(N,0,3)
  pre_censoring=pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
  tcens=(pre_censoring<pre_time) # censoring indicator
  delta=1-tcens
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  delta = delta[order(time)]
  facility=facility[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  dat = list(delta = delta, time= time, facility = facility,z = z)
  return(dat)
}

# Newton-Raphson method
TimeVarying_NR <- function(time, delta, z, facility = NULL, knot = 5, M_stop = 500, rate = 0.001, tol = 10^(-6)){
  p = ncol(z)
  N = nrow(z)
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(knot-4))/(knot-3))
  bs7=splines::bs(time,df=knot, knot=knot_set, intercept=TRUE, degree=3)
  bs8=matrix(bs7, nrow=N) # for build beta_t
  if(is.null(facility)){
    facility = rep(1,N)
  }
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  theta_stratify=matrix( rep(0, knot*p), nrow=p)  #dim P*knot
  theta_NR=matrix(rep(constant_beta, knot), nrow=p, byrow=FALSE) # initialize theta
  likelihood_NR_all=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_NR)/N # initial likelihood
  key=0
  repeat{
    key=key+1
    temp=ddloglik(knot,facility,delta,z,bs8,theta_NR,length(unique(facility)))
    dist=matrix(solve(temp$GVG)%*%matrix(temp$GR_test,ncol=1), nrow=p, byrow=TRUE) # NR update
    ## to find a gamma s.t. llk > llk_old
    gamma=1
    theta_temp=theta_NR+gamma*dist
    likelihood_stratify=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_temp)/N
    while(likelihood_stratify<likelihood_NR_all[key]+rate*gamma){
      gamma=gamma/2 #gradually decrease step size
      theta_temp=theta_NR+gamma*dist
      likelihood_stratify=dloglik_likelihood_gradient(knot,facility,delta,z,bs8,theta_temp)/N
    }
    
    theta_NR=theta_NR+gamma*dist
    likelihood_stratify=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_NR)/N  
    likelihood_NR_all=c(likelihood_NR_all, likelihood_stratify)
    if (key>=2) {
      llk.diff = likelihood_NR_all[key]-likelihood_NR_all[key-1]
      if(abs(llk.diff) < tol |max(abs(dist))<tol) break
    }
    if(key>=M_stop) break
  }
  beta_NR=matrix(rep(0, N*p), nrow=N)
  for (j in 1:p) {
    beta_NR[,j]=bs8%*%theta_NR[j,]
  }
  
  test_NR=rep(0,p)
  constrain=-diff(diag(knot*1),differences=1)
  j=0
  repeat{
    j=j+1
    theta_NR_j= theta_NR[j,]
    L2=solve(temp$GVG)[((j-1)*knot+1):((j-1)*knot+knot),((j-1)*knot+1):((j-1)*knot+knot)]
    test_contrast=t(constrain%*%theta_NR_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_NR_j)
    test_NR[j]=1-pchisq(test_contrast, (knot-1))
    if(j==p) break
  }
  
  rslt =list(theta = theta_NR, beta = beta_NR, llk_all = likelihood_NR_all, llk = likelihood_stratify, pvalue = test_NR, basis = bs8)
  return(rslt)
}

# GDboost method
TimeVarying_GDboost <- function(time, delta, z, facility = NULL, knot= 5, M_stop = 500, rate = 0.01, tol = 10^(-6)){
  p = ncol(z)
  N = nrow(z)
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(knot-4))/(knot-3))
  bs7=splines::bs(time,df=knot, knot=knot_set, intercept=TRUE, degree=3)
  bs8=matrix(bs7, nrow=N) # for build beta_t
  if(is.null(facility)){
    facility = rep(1,N)
  }
  m_stratify=0
  theta_GD=matrix(rep(0, knot*p), nrow=p)  #dim P*knot
  likelihood_GD_all=dloglik_likelihood_gradient(knot,facility,delta,z,bs8,theta_GD)/N
  stage2_key=FALSE
  track = 5 
  while (stage2_key==FALSE){
    m_stratify=m_stratify+1
    result=GDboost_stratify(knot,rate,facility,delta,z,bs8,theta_GD)
    theta_GD=result$theta
    likelihood_stratify=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_GD)/N
    likelihood_GD_all=c(likelihood_GD_all, likelihood_stratify)
    if (m_stratify>=(10+track)) {
      llk.diff = likelihood_GD_all[m_stratify]-likelihood_GD_all[m_stratify-1]
      llk.diff2 = likelihood_GD_all[m_stratify]-likelihood_GD_all[1]
      if(abs(llk.diff/llk.diff2) < tol) {
        stage2_key=TRUE
        break
      }
    }
    if (m_stratify==M_stop){
      stage2_key=TRUE
      break
    }
  } #end while
  beta_GD = matrix(rep(0, N*p), nrow=N)
  for (j in 1:p) {
    beta_GD[,j]=bs8%*%theta_GD[j,]
  }
  
  temp = ddloglik(knot,facility,delta,z,bs8,theta_GD,length(unique(facility)))
  test_GD=rep(0,p)
  constrain=-diff(diag(knot*1),differences=1)
  j=0
  repeat{
    j=j+1
    theta_GD_j= theta_GD[j,]
    L2=solve(temp$GVG)[((j-1)*knot+1):((j-1)*knot+knot),((j-1)*knot+1):((j-1)*knot+knot)]
    test_contrast=t(constrain%*%theta_GD_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_GD_j)
    test_GD[j]=1-pchisq(test_contrast, (knot-1))
    if(j==p) break
  }
  rslt =list(theta = theta_GD, beta = beta_GD, llk_all = likelihood_GD_all, llk = likelihood_stratify, pvalue = test_GD, basis = bs8)
  return(rslt)
}
