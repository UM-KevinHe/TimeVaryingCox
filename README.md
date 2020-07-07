# TimeVaryingCox

## example.R
Include an example of how to use functions **TimeVarying_NR()** and **TimeVarying_GDboost()**.
Estimated running time is about 0.5 hour. 


## TimeVarying_NR(time, delta, z, facility, knot, M_stop, rate, tol)
## Method: Newton-Raphson (self-implemented)

### Arguments:
- time: survival time
- delta: event indicator
- z: covariates
- facility: stratification variable
- knot: number of B spline bases
- M_stop: maximum number of iterations
- rate: step size
- tol: convergence criterion

### Output:
A *list* object:
- theta: B spline coefficients
- beta: varying coeffcients
- llk_all: likelihood for each iteration
- llk: final likelihood
- pvalue: p values for contant effect tests
- basis: B spline basis

## TimeVarying_GDboost(time, delta, z, facility, knot, M_stop, rate, track, tol)
### Method: GDBoosting 

### Arguments:
- time: survival time
- delta: event indicator
- z: covariates
- facility: stratification variable
- knot: number of B spline bases
- M_stop: maximum number of iterations
- rate: step size
- tol: convergence criterion

### Output:
A *list* object:
- theta: B spline coefficients
- beta: varying coeffcients
- llk_all: likelihood for each iteration
- llk: final likelihood
- pvalue: p values for contant effect tests
- basis: B spline basis

## R package BSATV
**BSATV** includes an R package for this project. 

