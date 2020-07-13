library(survival)

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


# generate covariates with different cov structure
covariate_sim <- function(n, p, z.type = "Ind", z.para) {
  z = matrix(NA, n, p)
  if (z.type == "Ind") {
    for (i in 1:p) {
      z[, i] = runif(n, 0, 1)
    }
  }
  if (z.type == "AR1") {
    z_mat = diag(p)
    z_mat = z.para^abs(row(z_mat) - col(z_mat))
    z = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = z_mat)
  }
  if (z.type == "CS") {
    z_mat = matrix(z.para, p, p)
    diag(z_mat) = rep(1, p)
    z = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = z_mat)
  }
  colnames(z) = paste("z", 1:p, sep = "")
  return(z)
}


#' Sampling within stratification.
#'
#' @param size A vector of sizes of each stratum.
#' @param prop Proportion of sample size to original size.
#' @return A vector of sample indices.
#' @export
#' @examples
#' stra_sampling(c(3,5,3), prop = 1)
# stratified sampling
stra_sampling <- function(size, prop) {
  m <- length(size)
  mark <- cumsum(size)
  total <- sum(size)
  sample.id <- NULL
  sample.id <- c(sample.id, sample.int(size[1], size = floor(size[1]/prop), replace = FALSE))
  for (i in 2:m) {
    sample.id <- c(sample.id, mark[i - 1] + sample.int(size[i], size = floor(size[i]/prop), replace = FALSE))
  }
  return(sample.id)
}

#' Simulate survival data from the stratified Cox model with time-varying coefficients.
#'
#' @param nF Number of strata.
#' @param alph Average size of strata.
#' @param seed Random seed.
#' @param censor.max Maximum value in the censoring time distribution.
#' @param end.time Censoring time.
#' @param fun.list A list of coefficient functions.
#' @param z.type Covariance type: 'Ind' (independent), 'AR(1)' (AR(1)), and 'CS' (compound symmetry)
#' @param z.para Correlation coefficient for covariance matrix.
#' @return A list containing the stratified survival data.
#' \describe{
#'   \item{delta}{event indicator}
#'   \item{time}{survival time}
#'   \item{facility}{stratum ID}
#'   \item{z}{covariates}
#' }
#'
#' @export
#' @examples
#' z.type = 'Ind' # 'Ind', 'AR1', 'CS'
#' z.para = 0.5
#' f1<-function(x){
#'  (-x^2+3)
#' }
#' f1 =  Vectorize(f1)
#' f2<-function(x){
#'  2*log(x+0.01)
#'  }
#' f2 = Vectorize(f2)
#' f3<-function(x){
#'   2*(-3/(x+1)+1)
#' }
#' f3 = Vectorize(f3)
#' fun.list=list(f1,f2,f3)
#' dat = sim_data_fun(nF = 50, alph=50,seed=12,censor.max=30,end.time=3,fun.list=fun.list, z.type = z.type, z.para = z.para)
# data simulation simulate multi-variable functions by function list
sim_data_fun <- function(nF = 50, alph = 50, seed = 12, censor.max = 30, end.time = 3, fun.list = fun.list, z.type = z.type,
                         z.para = z.para) {
  set.seed(seed)
  p = length(fun.list)
  n_f = rpois(nF, lambda = alph)  #sample size: a lot of 00
  N = sum(n_f)
  ############ generate data########################
  gamma = rnorm(nF, mean = 0, sd = 0.5)
  gamma_subject = rep(gamma, n_f)
  z = covariate_sim(N, p, z.type = z.type, z.para = z.para)
  U = runif(N, 0, 1)
  F_pre = 1:nF
  facility = rep(F_pre, n_f)
  time = rep(0, N)
  kerf <- function(x, i) {
    ker = 0
    for (j in 1:p) {
      ker = ker + fun.list[[j]](x) * z[i, j]
    }
    ker
  }
  for (i in 1:N) {
    f = function(t) {
      integrand <- function(x) {
        gamma0 * exp(gamma_subject[i] + kerf(x, i))
      }  #3*sin(3*pi*x/4)*(x<3)*z[i,1]
      alph = integrate(integrand, lower = 0, upper = t)$value
      alph + log(1 - U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f, lower = 0, upper = 4), silent = TRUE))
    if (class(r1) == "try-error") {
      time[i] = 4
    } else time[i] = uniroot(f, lower = 0, upper = 4)$root
  }
  censoring = runif(N, 1, censor.max)
  censoring = censoring * (censoring < 3) + 3 * (censoring >= 3)
  tcens = (censoring < time)  # censoring indicator
  delta = 1 - tcens
  time = time * (delta == 1) + censoring * (delta == 0)
  # mean(delta) # 0.8
  delta = delta[order(time)]
  z = z[order(time), ]
  facility = facility[order(time)]
  time = time[order(time)]
  return(list(delta = delta, z = z, time = time, facility = facility))
}


#' Simulate survival data from the stratified Cox model with time-varying coefficients.
#'
#' @param nF Number of strata.
#' @param p Number of covariates.
#' @param seed Random seed.
#' @return A list containing the stratified survival data.
#' \describe{
#'   \item{delta}{event indicator}
#'   \item{time}{survival time}
#'   \item{facility}{stratum ID}
#'   \item{z}{covariates}
#' }
#' @export
#' @examples
#' dat1 = sim_data(nF = 20, p = 5, seed = 1)
# data simulation
sim_data <- function(nF = 20, p = 5, seed = 1) {
  set.seed(seed)
  n_f = rpois(nF, lambda = 500)  #sample size for each facility
  N = sum(n_f)
  gamma = rnorm(nF, mean = 0, sd = 0.5)
  range(gamma)
  gamma_subject = rep(gamma, n_f)
  F_pre = 1:nF
  facility = rep(F_pre, n_f)
  ############ generate data########################
  Sigma_z1 <- AR1(0.6, p)
  z = mvtnorm::rmvnorm(N, mean = rep(0, p), sigma = Sigma_z1)
  z_012_rare = function(x) {
    # U=runif(2, 0.97,1)
    U = runif(2, 0.9, 0.975)
    U = U[order(U)]
    x2 = quantile(x, prob = U)
    x3 = x
    x3[x < x2[1]] = 0
    x3[x > x2[2]] = 2
    x_index = (((x < x2[1]) + (x > x2[2])) < 1)
    x3[x_index] = 1
    return(x3)
  }
  z[, 1:3] = apply(z[, 1:3], 2, z_012_rare)
  MAF = apply(z, 2, mean)
  range(MAF)  #0.0769688 0.1778108
  U = runif(N, 0, 1)
  pre_time = rep(0, N)
  for (i in 1:(N)) {
    f = function(t) {
      integrand <- function(x) {
        0.5 * exp(gamma_subject[i] + z[i, 1] - z[i, 3] + sin(3 * pi * x/4) * (x < 3) * z[i, 2] - (x/3)^2 * exp(x/2) *
                    (x < 3) * z[i, 4] + z[i, 5])
      }
      Lambda = integrate(integrand, lower = 0, upper = t)$value
      Lambda + log(1 - U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f, lower = 0, upper = 4), silent = TRUE))
    if (class(r1) == "try-error") {
      pre_time[i] = 4
    } else pre_time[i] = uniroot(f, lower = 0, upper = 4)$root
  }
  
  pre_censoring = runif(N, 0, 3)
  pre_censoring = pre_censoring * (pre_censoring < 3) + 3 * (pre_censoring >= 3)
  tcens = (pre_censoring < pre_time)  # censoring indicator
  delta = 1 - tcens
  time = pre_time * (delta == 1) + pre_censoring * (delta == 0)
  delta = delta[order(time)]
  facility = facility[order(time)]
  z = z[order(time), ]
  time = time[order(time)]
  dat = list(delta = delta, time = time, facility = facility, z = z)
  return(dat)
}



#' Fit stratified time-varying Cox model using Newton-Raphson algorithm.
#'
#' @param time Survival time.
#' @param delta Event indicator.
#' @param z Covariates.
#' @param facility Stratum ID.
#' @param knot Number of B spline bases.
#' @param M_stop Maximum number of iteration steps.
#' @param rate Step size.
#' @param tol Convergen criterion.
#' @return A list of fitting results.
#' \describe{
#'    \item{theta}{Estimated B spline coefficients}
#'    \item{beta}{Estimated time-varying coefficient functions}
#'    \item{llk_all}{Fitted partial likelihood in each step}
#'    \item{likelihood_stratify}{likelihood_stratify}
#'    \item{pvalue}{P values for testing constant coefficients}
#'    \item{basis}{B spline bases}
#' }
#' @export
#' @examples
#' \dontrun{
#' delta = dat1$delta
#' facility = dat1$facility
#' z = dat1$z
#' time = dat1$time
#' #fit_NR = TimeVarying_NR(time, delta, z, facility, knot = 6, M_stop = 500, rate = 0.001, tol = 10^(-6))
#' }
# Newton-Raphson method
TimeVarying_NR <- function(time, delta, z, facility = NULL, knot = 5, M_stop = 500, rate = 0.001, tol = 10^(-6)) {
  p = ncol(z)
  N = nrow(z)
  if (is.null(facility)) {
    facility = rep(1, N)
  }
  delta = delta[order(time)]
  facility = facility[order(time)]
  z = z[order(time), ]
  time = time[order(time)]
  time2 = time[delta == 1]
  knot_set = quantile(time2, prob = seq(1:(knot - 4))/(knot - 3))
  bs7 = splines::bs(time, df = knot, knot = knot_set, intercept = TRUE, degree = 3)
  bs8 = matrix(bs7, nrow = N)  # for build beta_t
  
  vfit <- survival::coxph(Surv(time, delta) ~ z)
  constant_beta <- vfit$coef
  theta_stratify = matrix(rep(0, knot * p), nrow = p)  #dim P*knot
  theta_NR = matrix(rep(constant_beta, knot), nrow = p, byrow = FALSE)  # initialize theta
  likelihood_NR_all = dloglik_likelihood_stratify(knot, facility, delta, z, bs8, theta_NR)/N  # initial likelihood
  key = 0
  repeat {
    key = key + 1
    temp = ddloglik(knot, facility, delta, z, bs8, theta_NR, length(unique(facility)))
    dist = matrix(solve(temp$GVG) %*% matrix(temp$GR_test, ncol = 1), nrow = p, byrow = TRUE)  # NR update
    ## to find a gamma s.t. llk > llk_old
    gamma = 1
    theta_temp = theta_NR + gamma * dist
    likelihood_stratify = dloglik_likelihood_stratify(knot, facility, delta, z, bs8, theta_temp)/N
    while (likelihood_stratify < likelihood_NR_all[key] + rate * gamma) {
      gamma = gamma/2  #gradually decrease step size
      theta_temp = theta_NR + gamma * dist
      likelihood_stratify = dloglik_likelihood_gradient(knot, facility, delta, z, bs8, theta_temp)/N
    }
    
    theta_NR = theta_NR + gamma * dist
    likelihood_stratify = dloglik_likelihood_stratify(knot, facility, delta, z, bs8, theta_NR)/N
    likelihood_NR_all = c(likelihood_NR_all, likelihood_stratify)
    if (key >= 2) {
      llk.diff = likelihood_NR_all[key] - likelihood_NR_all[key - 1]
      if (abs(llk.diff) < tol | max(abs(dist)) < tol)
        break
    }
    if (key >= M_stop)
      break
  }
  beta_NR = matrix(rep(0, N * p), nrow = N)
  for (j in 1:p) {
    beta_NR[, j] = bs8 %*% theta_NR[j, ]
  }
  
  test_NR = rep(0, p)
  constrain = -diff(diag(knot * 1), differences = 1)
  j = 0
  repeat {
    j = j + 1
    theta_NR_j = theta_NR[j, ]
    L2 = solve(temp$GVG)[((j - 1) * knot + 1):((j - 1) * knot + knot), ((j - 1) * knot + 1):((j - 1) * knot + knot)]
    test_contrast = t(constrain %*% theta_NR_j) %*% solve(constrain %*% L2 %*% t(constrain)) %*% (constrain %*%
                                                                                                    theta_NR_j)
    test_NR[j] = 1 - pchisq(test_contrast, (knot - 1))
    if (j == p)
      break
  }
  
  rslt = list(theta = theta_NR, beta = beta_NR, llk_all = likelihood_NR_all, llk = likelihood_stratify, pvalue = test_NR,
              basis = bs8)
  return(rslt)
}



#' Fit stratified time-varying Cox model using block-wise steepest ascent procedure.
#'
#' @param time Survival time.
#' @param delta Event indicator.
#' @param z Covariates.
#' @param facility Stratum ID.
#' @param knot Number of B spline bases.
#' @param M_stop Maximum number of iteration steps.
#' @param rate Step size.
#' @param tol Convergen criterion.
#' @return A list of fitting results.
#' \describe{
#'    \item{theta}{Estimated B spline coefficients}
#'    \item{beta}{Estimated time-varying coefficient functions}
#'    \item{llk_all}{Fitted partial likelihood in each step}
#'    \item{likelihood_stratify}{likelihood_stratify}
#'    \item{pvalue}{P values for testing constant coefficients}
#'    \item{basis}{B spline bases}
#' }
#' @export
#' @examples
#' \dontrun{
#' delta = dat1$delta
#' facility = dat1$facility
#' z = dat1$z
#' time = dat1$time
#' fit_GDboost = TimeVarying_GDboost(time, delta, z, facility, knot = 6, M_stop = 500, rate = 0.001, tol = 10^(-6))
#' }
# GDboost method
TimeVarying_GDboost <- function(time, delta, z, facility = NULL, knot = 5, M_stop = 500, rate = 0.01, tol = 10^(-6)) {
  p = ncol(z)
  N = nrow(z)
  if (is.null(facility)) {
    facility = rep(1, N)
  }
  delta = delta[order(time)]
  facility = facility[order(time)]
  z = z[order(time), ]
  time = time[order(time)]
  time2 = time[delta == 1]
  knot_set = quantile(time2, prob = seq(1:(knot - 4))/(knot - 3))
  bs7 = splines::bs(time, df = knot, knot = knot_set, intercept = TRUE, degree = 3)
  bs8 = matrix(bs7, nrow = N)  # for build beta_t
  
  m_stratify = 0
  theta_GD = matrix(rep(0, knot * p), nrow = p)  #dim P*knot
  likelihood_GD_all = dloglik_likelihood_gradient(knot, facility, delta, z, bs8, theta_GD)/N
  stage2_key = FALSE
  track = 5
  while (stage2_key == FALSE) {
    m_stratify = m_stratify + 1
    result = GDboost_stratify(knot, rate, facility, delta, z, bs8, theta_GD)
    theta_GD = result$theta
    likelihood_stratify = dloglik_likelihood_stratify(knot, facility, delta, z, bs8, theta_GD)/N
    likelihood_GD_all = c(likelihood_GD_all, likelihood_stratify)
    if (m_stratify >= (10 + track)) {
      llk.diff = likelihood_GD_all[m_stratify] - likelihood_GD_all[m_stratify - 1]
      llk.diff2 = likelihood_GD_all[m_stratify] - likelihood_GD_all[1]
      if (abs(llk.diff/llk.diff2) < tol) {
        stage2_key = TRUE
        break
      }
    }
    if (m_stratify == M_stop) {
      stage2_key = TRUE
      break
    }
  }  #end while
  beta_GD = matrix(rep(0, N * p), nrow = N)
  for (j in 1:p) {
    beta_GD[, j] = bs8 %*% theta_GD[j, ]
  }
  
  temp = ddloglik(knot, facility, delta, z, bs8, theta_GD, length(unique(facility)))
  test_GD = rep(0, p)
  constrain = -diff(diag(knot * 1), differences = 1)
  j = 0
  repeat {
    j = j + 1
    theta_GD_j = theta_GD[j, ]
    L2 = solve(temp$GVG)[((j - 1) * knot + 1):((j - 1) * knot + knot), ((j - 1) * knot + 1):((j - 1) * knot + knot)]
    test_contrast = t(constrain %*% theta_GD_j) %*% solve(constrain %*% L2 %*% t(constrain)) %*% (constrain %*%
                                                                                                    theta_GD_j)
    test_GD[j] = 1 - pchisq(test_contrast, (knot - 1))
    if (j == p)
      break
  }
  rslt = list(theta = theta_GD, beta = beta_GD, llk_all = likelihood_GD_all, llk = likelihood_stratify, pvalue = test_GD,
              basis = bs8)
  return(rslt)
}


TimeVarying_SGD <- function(time, delta, z, facility = NULL, knot = 5, M_stop = 500, rate = 0.01, tol = 10^(-6), seed = 1) {
  p = ncol(z)
  N = nrow(z)
  if (is.null(facility)) {
    facility = rep(1, N)
  }
  delta = delta[order(time)]
  facility = facility[order(time)]
  z = z[order(time), ]
  time = time[order(time)]
  time2 = time[delta == 1]
  knot_set = quantile(time2, prob = seq(1:(knot - 4))/(knot - 3))
  bs7 = splines::bs(time, df = knot, knot = knot_set, intercept = TRUE, degree = 3)
  bs8 = matrix(bs7, nrow = N)  # for build beta_t
  m = as.vector(sapply(split(facility,factor(facility)),length))
  m_stratify=0
  theta_stratify=matrix( rep(0, knot*p), nrow=p)  #dim P*knot
  
  G=rep(0,p)
  gamma=0.01 #0.01 performs better (better than 1, 0.1, 0.001, 0.0001), and is recommended by https://arxiv.org/pdf/1609.04747.pdf
  stage2_key = FALSE
  likelihood_stochastic_gradient_all=NULL
  while (stage2_key==FALSE){
    m_stratify=m_stratify+1
    set.seed(seed*1000+m_stratify)
    bootsample3 <- stra_sampling(m,10)
    bootsample3 = sort(bootsample3)
    time_sub4=time[bootsample3]
    delta_sub4=delta[bootsample3]
    facility_sub4=facility[bootsample3]
    b_spline_sub4=bs8[bootsample3,]
    z_sub4=z[bootsample3,]
    
    result=GDboost_gradient(knot,rate,facility_sub4,delta_sub4,z_sub4,b_spline_sub4,theta_stratify)
    if (m_stratify>1){
      theta_stratify=theta_stratify-diag(gamma/sqrt(G))%*%matrix(result$L1,nrow=p,byrow = TRUE)
    }
    if (m_stratify==1){
      theta_stratify=theta_stratify-gamma*matrix(result$L1,nrow=p,byrow = TRUE)
    }
    G_temp=rep(0,p)
    for (j in 1:p){
      G_temp[j]=sum((result$L1[((j-1)*knot+1):(j*knot)])**2)/knot
    }
    G=G+G_temp
    likelihood_stratify=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_stratify)/N  ###a function written in function.R
    likelihood_stochastic_gradient_all=c(likelihood_stochastic_gradient_all, likelihood_stratify)
    beta_stochastic_gradient=matrix(rep(0, N*p), nrow=N)
    diff_stochastic_gradient=0
    for (j in 1:p) {
      beta_stochastic_gradient[,j]=bs8%*%theta_stratify[j,]
    }
    track = 5
    if (m_stratify>=(10+track)) {
      llk.diff.all=NULL
      for (key in 1:track) {
        llk.diff.all = c(llk.diff.all, likelihood_stochastic_gradient_all[m_stratify-key+1]-likelihood_stochastic_gradient_all[m_stratify-key])
      }
      llk.diff=max(llk.diff.all)
      llk.diff_null = likelihood_stochastic_gradient_all[m_stratify]-likelihood_stochastic_gradient_all[1]
      if(abs(llk.diff/llk.diff_null) < tol) {
        stage2_key=TRUE
        break
      }
    }
    if (m_stratify==M_stop){
      stage2_key=TRUE
      break
    }
  } #end while
  theta_SGD = theta_stratify
  temp = ddloglik(knot, facility, delta, z, bs8, theta_SGD, length(unique(facility)))
  beta_SGD = matrix(rep(0, N * p), nrow = N)
  for (j in 1:p) {
    beta_SGD[, j] = bs8 %*% theta_SGD[j, ]
  }
  test_SGD = rep(0, p)
  constrain = -diff(diag(knot * 1), differences = 1)
  j = 0
  repeat {
    j = j + 1
    theta_SGD_j = theta_SGD[j, ]
    L2 = solve(temp$GVG)[((j - 1) * knot + 1):((j - 1) * knot + knot), ((j - 1) * knot + 1):((j - 1) * knot + knot)]
    test_contrast = t(constrain %*% theta_SGD_j) %*% solve(constrain %*% L2 %*% t(constrain)) %*% (constrain %*%theta_SGD_j)
    test_SGD[j] = 1 - pchisq(test_contrast, (knot - 1))
    if (j == p)
      break
  }
  rslt = list(theta = theta_SGD, beta = beta_SGD, llk_all = likelihood_stochastic_gradient_all, llk = likelihood_stratify, pvalue = test_SGD,
              basis = bs8)
  return(rslt)
}

