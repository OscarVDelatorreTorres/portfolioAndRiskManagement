
# portfolio parameters functions:
# This R script includes all the parameter functions for portfolio optimization in the fPortfolio package.
# More specifically, it includes functions for the expected returns and covariance matrix of the assets used
# in the setEstimator() fPortfolio function.
# setEstimator(portfoliospec2)=myestimator

# Version control:
# V 1.0: Oscar e la Torre-Torres (original script)

# Required libraries upload or instalation:
if (!require(fPortfolio)) {install.packages('fPortfolio')
  library(fPortfolio)} else {library(fPortfolio)}
if (!require(rmgarch)) {install.packages('rmgarch')
  library(rmgarch)} else {library(rmgarch)}
if (!require(tseries)) {install.packages('tseries')
  library(tseries)} else {library(tseries)}
if (!require(xts)) {install.packages('xts')
  library(xts)} else {library(xts)}
# Time-fixed:====

timeFixedCovariance=function(Returns,portfoliospec){
  rendimientos=as.data.frame(Returns)
  
  covarianza=as.matrix(cov(rendimientos))
  medias=as.matrix(colMeans(rendimientos))
  
  list(mu=medias,
       Sigma=covarianza)
}

# A: symmetric DCC models====
# A1 DCC symetric GARCH models====

# DCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian GARCH:====

gaussianDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)

  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm")
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GARCH:====

gaussianDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and GED GARCH:====

gaussianDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian GARCH:====

gaussianDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GARCH:====

gaussianDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED GARCH:====

gaussianDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-T correlation and Gaussian GARCH:====

stdDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent GARCH:====

stdDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and GED GARCH:====

stdDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-t correlation and sGaussian GARCH:====

gaussianDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent GARCH:====

stdDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")

  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and skewed GED GARCH:====

stdDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# DCCGARCH-(0,0)ARMA La Place correlation and Gaussian GARCH:====

laplaceDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")

  print("Estimating La Place correlation and Gaussian ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent GARCH:====

laplaceDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and GED GARCH:====

laplaceDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA La Place correlation and sGaussian GARCH:====

laplaceDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent GARCH:====

laplaceDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and skewed GED GARCH:====

laplaceDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed GED ARMA(0,0) GARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}

# A2 DCC EGARCH models====

# DCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian EGARCH:====

gaussianDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent EGARCH:====

gaussianDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and GED EGARCH:====

gaussianDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian EGARCH:====

gaussianDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent EGARCH:====

gaussianDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED EGARCH:====

gaussianDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-T correlation and Gaussian EGARCH:====

stdDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent EGARCH:====

stdDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and GED EGARCH:====

stdDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-t correlation and sGaussian EGARCH:====

gaussianDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent EGARCH:====

stdDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and skewed GED EGARCH:====

stdDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# DCCGARCH-(0,0)ARMA La Place correlation and Gaussian EGARCH:====

laplaceDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and Gaussian ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent EGARCH:====

laplaceDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and GED EGARCH:====

laplaceDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA La Place correlation and sGaussian EGARCH:====

laplaceDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent EGARCH:====

laplaceDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and skewed GED EGARCH:====

laplaceDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed GED ARMA(0,0) EGARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}

# A3 DCC GJRGARCH models====

# DCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian GJRGARCH:====

gaussianDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GJRGARCH:====

gaussianDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) gjrGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and GED GJRGARCH:====

gaussianDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian GJRGARCH:====

gaussianDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GJRGARCH:====

gaussianDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED GJRGARCH:====

gaussianDccSGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="DCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-T correlation and Gaussian GJRGARCH:====

stdDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent GJRGARCH:====

stdDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and GED GJRGARCH:====

stdDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA Student-t correlation and sGaussian GJRGARCH:====

gaussianDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and sStudent GJRGARCH:====

stdDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA Student-t correlation and skewed GED GJRGARCH:====

stdDccSGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="DCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# DCCGARCH-(0,0)ARMA La Place correlation and Gaussian GJGARCH:====

laplaceDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and Gaussian ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent GJRGARCH:====

laplaceDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and GED GJRGARCH:====

laplaceDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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


# DCCGARCH-(0,0)ARMA La Place correlation and sGaussian GJRGARCH:====

laplaceDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and sStudent GJGARCH:====

laplaceDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
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

# DCCGARCH-(0,0)ARMA La Place correlation and skewed GED GJRGARCH:====

laplaceDccSGEDGGJRARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="DCC")
  print("Estimating La Place correlation and skewed GED ARMA(0,0) GJRGARCH DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# B: asymmetric DCC models====
# B1 aDCC symetric GARCH models====

# aDCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian GARCH:====

gaussianaDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm")
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GARCH:====

gaussianaDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and GED GARCH:====

gaussianaDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) GARCH asymmetricDCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian GARCH:====

gaussianaDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GARCH:====

gaussianaDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED GARCH:====

gaussianaDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-T correlation and Gaussian GARCH:====

stdaDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent GARCH:====

stdaDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and GED GARCH:====

stdaDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-t correlation and sGaussian GARCH:====

gaussianaDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent GARCH:====

stdaDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and skewed GED GARCH:====

stdaDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# aDCCGARCH-(0,0)ARMA La Place correlation and Gaussian GARCH:====

laplaceaDccGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Gaussian ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent GARCH:====

laplaceaDccTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and GED GARCH:====

laplaceaDccGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and GED ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA La Place correlation and sGaussian GARCH:====

laplaceaDccSGaussianGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) GARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent GARCH:====

laplaceaDccSTStudentGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) GARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and skewed GED GARCH:====

laplaceaDccSGEDGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed GED ARMA(0,0) GARCH DCC-GARCH asymmetric covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}

# B2 aDCC EGARCH models====

# aDCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian EGARCH:====

gaussianaDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) EGARCH DCC-GARCH asymmetric covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent EGARCH:====

gaussianaDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and GED EGARCH:====

gaussianaDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian EGARCH:====

gaussianaDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent EGARCH:====

gaussianaDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED EGARCH:====

gaussianaDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-T correlation and Gaussian EGARCH:====

stdaDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent EGARCH:====

stdaDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and GED EGARCH:====

stdaDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-t correlation and sGaussian EGARCH:====

gaussianaDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent EGARCH:====

stdaDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and skewed GED EGARCH:====

stdaDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# aDCCGARCH-(0,0)ARMA La Place correlation and Gaussian EGARCH:====

laplaceaDccGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Gaussian ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent EGARCH:====

laplaceaDccTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and GED EGARCH:====

laplaceaDccGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA La Place correlation and sGaussian EGARCH:====

laplaceaDccSGaussianEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) EGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent EGARCH:====

laplaceaDccSTStudentEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and skewed GED EGARCH:====

laplaceaDccSGEDEGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="eGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed GED ARMA(0,0) EGARCH asymmetric DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}

# B3 aDCC GJRGARCH models====

# aDCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian GJRGARCH:====

gaussianaDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Gaussian ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GJRGARCH:====

gaussianaDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and Student-t ARMA(0,0) gjrGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and GED GJRGARCH:====

gaussianaDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Gaussian correlation and sGaussian GJRGARCH:====

gaussianaDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and sStudent GJRGARCH:====

gaussianaDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed Student-t ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Gaussian correlation and skewed GED GJRGARCH:====

gaussianaDccSGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvnorm",model="aDCC")
  print("Estimating Gaussian correlation and skewed GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-T correlation and Gaussian GJRGARCH:====

stdaDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Gaussian ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent GJRGARCH:====

stdaDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and Student-t ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and GED GJRGARCH:====

stdaDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA Student-t correlation and sGaussian GJRGARCH:====

gaussianaDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and sStudent GJRGARCH:====

stdaDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed Student-t ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA Student-t correlation and skewed GED GJRGARCH:====

stdaDccSGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvt",model="aDCC")
  print("Estimating Student-t correlation and skewed GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}


# aDCCGARCH-(0,0)ARMA La Place correlation and Gaussian GJGARCH:====

laplaceaDccGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="norm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Gaussian ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent GJRGARCH:====

laplaceaDccTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="std"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and Student-t ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and GED GJRGARCH:====

laplaceaDccGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="ged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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


# aDCCGARCH-(0,0)ARMA La Place correlation and sGaussian GJRGARCH:====

laplaceaDccSGaussianGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="snorm"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Gaussian ARMA(0,0) GJRGARCH with Gaussian asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and sStudent GJGARCH:====

laplaceaDccSTStudentGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sstd"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")
  print("Estimating La Place correlation and skewed Student-t ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
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

# aDCCGARCH-(0,0)ARMA La Place correlation and skewed GED GJRGARCH:====

laplaceaDccSGEDGJRGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model="gjrGARCH", garchOrder=c(1, 1)),
    mean.model=list(armaOrder=c(0, 0), include.mean=TRUE),
    distribution.model="sged"
  )
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                      distribution = "mvlaplace",model="aDCC")

  print("Estimating La Place correlation and skewed GED ARMA(0,0) GJRGARCH asymmetric DCC-GARCH covariance matrix...")
  # fit the DCC spec to returns
  dcc_fit <- dccfit(dcc_spec, rendimientos)
  
  forecastGARCH=dccforecast(dcc_fit, n.ahead=1)
  
  # 1 month forecast (i.e 4 weeks ahead)
  covarianza <- rcov(forecastGARCH)[[1]]
  covarianza=covarianza[,,1]
  # mean vector:
  medias=as.matrix(forecastGARCH@mforecast[["mu"]])
  #---
  # Output list:
  list(mu=medias,
       Sigma=covarianza)
}

