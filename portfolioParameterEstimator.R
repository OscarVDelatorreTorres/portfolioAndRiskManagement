
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

# Input parameters:

#garchModel: univariate model GARCH model
#garchLLF: Univariate GARCH model LLF.
#armaOrder: The order of the ARMA model.
#garchOrder: the order of the GARCH model.
#varDCC: if FALSE, the DCC model is estimated on the residuals. If TRUE, 
# the DCC model is estimated from a VAR model's residuals.
#varLag: The maximum lag of the VAR model
#dccPdf: the multivariate pdf of the DCC correlation estimation.
#dccModel: the DCC model to esitmate. Either "dccModel" for a symmetric models
# or "adccModel" forn an assymetric one.

dccGARCH=function(Returns,portfoliospec){
  
  rendimientos=as.data.frame(Returns)

  #---
  # GARCH(1,1) Specification
  garch_spec <- ugarchspec(
    variance.model=list(model=garchModel, garchOrder=garchOrder),
    mean.model=list(armaOrder=armaOrder, include.mean=TRUE),
    distribution.model=garchLLF)
  
  # create multispec--a set of GARCH(1,1) specifications on each series
  ms <- multispec(replicate(ncol(rendimientos), garch_spec))
  
  # turn multispec into a DCC spec
  
  if (varOrder>0){
  
    dcc_spec <- dccspec(ms, VAR = TRUE, robust = TRUE, lag = 1, lag.max = NULL,
                        distribution = dccPdf,model=dccModel)

    dccCorrelModel=paste0("Estimating ",dccPdf," ",dccModel," correlation with",
                          " VARMA (",varOrder,",",varOrder,")-",
                          garchLLF," ",garchModel,
                          " (",paste(as.character(garchOrder),collapse=","),") ",
                          "covariance matrix...")   
    
  } else {
 
    dcc_spec <- dccspec(ms, VAR = FALSE, robust = TRUE, lag = 1, lag.max = NULL,
                        distribution = dccPdf,model=dccModel)
    
    dccCorrelModel=paste0("Estimating ",dccPdf," ",dccModel," correlation with",
                          " ARMA (",paste(as.character(armaOrder),collapse=","),
                          ")-",
                          garchLLF," ",garchModel,
                          " (",paste(as.character(garchOrder),collapse=","),") ",
                          "covariance matrix...")    
    
  }
  
  
  print(dccCorrelModel)
  
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
       Sigma=covarianza,
       correlModel=dccCorrelModel)
}

