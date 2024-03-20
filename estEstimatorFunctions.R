
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

# DCCGARCH-(0,0)ARMA Gaussian correlation and Gaussian GARCH:====

gaussianDccGaussianGARCH=function(Returns,portfoliospec){
  
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
                      distribution = "mvnorm")
  
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
                      distribution = "mvnorm")
  
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
                      distribution = "mvnorm")
  
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
                      distribution = "mvnorm")
  
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
                      distribution = "mvnorm")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvt")
  
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
                      distribution = "mvlaplace")
  
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
                      distribution = "mvlaplace")
  
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
                      distribution = "mvlaplace")
  
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
                      distribution = "mvlaplace")
  
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
                      distribution = "mvlaplace")
  
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
                      distribution = "mvt")
  
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




#===============================================================================
# GOGARCH-(0,0) ARMA with Gaussian symmetric GARCH:====

GOGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)
  
  #---
  # GOGARCH(1,1) Specification

  goGarchSpec=gogarchspec(mean.model = list(model = "constant", robust = TRUE, 
                                            lag = 1, lag.max = NULL),
                          variance.model=list(model="sGARCH", 
                                              garchOrder=c(1, 1)),
                          distribution = "mvnorm"
                          )
                          
  
  print("GOGARCH ARMA(0,0) covariance matrix...")
  
  # fit the DCC spec to returns
  gogarch_fit <- gogarchfit(goGarchSpec, rendimientos)
  
  forecastGARCH=gogarchforecast(gogarch_fit, n.ahead=1)
  
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