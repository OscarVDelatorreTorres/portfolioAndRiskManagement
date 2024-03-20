
# portfolio parameters functions:
# This R script includes all the parameter functions for portfolio optimization in the fPortfolio package.
# More specifically, it includes functions for the expected returns and covariance matrix of the assets used
# in the setEstimator() fPortfolio function.
# setEstimator(portfoliospec2)=myestimator

if (!require(fPortfolio)) {install.packages('fPortfolio')
  library(fPortfolio)} else {library(fPortfolio)}
if (!require(rmgarch)) {install.packages('rmgarch')
  library(rmgarch)} else {library(rmgarch)}
if (!require(tseries)) {install.packages('tseries')
  library(tseries)} else {library(tseries)}
if (!require(xts)) {install.packages('xts')
  library(xts)} else {library(xts)}

# DCC-GARCH-(0,0) ARMA Gaussian===============================================

dccGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)

  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", submodel = "GARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm")
  
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Outpu list:
  list(mu=medias,
       Sigma=covarianza)
}

# DCC-GARCH-(1,1) ARMA===============================================

dccGARCHARMA=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)

  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", submodel = "GARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(1, 1), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm")
  
  print("Estimating Gaussian correlation and Gaussian ARMA(1,1) GARCH DCC-GARCH covariance matrix...")
  
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Outpu list:
  list(mu=medias,
       Sigma=covarianza)
}

# DCC-GARCH-(1,1) VARMA===============================================

dccGARCHVARMA=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)

  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", submodel = "GARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(1, 1), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = TRUE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm")
  print("Estimating Gaussian correlation and Gaussian VARMA(1,1) GARCH DCC-GARCH covariance matrix...")
  
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Outpu list:
  list(mu=medias,
       Sigma=covarianza)
}