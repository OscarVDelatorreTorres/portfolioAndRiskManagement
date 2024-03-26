# This is the min tracking error portfolio function for the ICTI 2023-2024 research proyect of
# Oscar De la Torre-Torres <https://www.oscardelatorretorres.com/icti2023>
# This routine estimates the minimum tracking error portfolio for a given benchmark portfolion
# In this case non-commodities such as avocado or green lemon.
# The idea of this function is to estimate the optimal portfolio results, along with the
# parameters of the minimum variance portfolio (minimum tracking error).
# The core idea is to estimate, given a Returns data matrix and a preselected
# portfolio parameters function, the theoretical portfolios and store the results in a database.

dccGARCHMinTrackingzErrorPortfolio=function(Returns,garchPortParams){
  
  armaOrder=garchPortParams$armaOrder
  garchModel=garchPortParams$garchModel
  garchLLF=garchPortParams$garchLLF
  garchOrder=garchPortParams$garchOrder
  varOrder=garchPortParams$varOrder
  dccPdf=garchPortParams$dccPdf
  dccModel=garchPortParams$dccModel
  
  # Extracting the asset names from the Return object:
  assetNames=colnames(Returns)[2:ncol(Returns)]
  
  # Portfolio frontiers======  
  # Conventional portfolio spec:
  portfoliospec=portfolioSpec(portfolio= list(
    weights = NULL, targetReturn = NULL,
    targetRisk = NULL, riskFreeRate = 0, nFrontierPoints = 100,
    status = NA))
  
  setEstimator(portfoliospec)=dccGARCH
  
  # Efficient frontiers estimation:
  print(paste0("Estimating efficient frontier for portfolio: ",experiment))
  portfolio1=portfolioFrontier(data=as.timeSeries(Returns),spec=portfoliospec)
  
  # Output efficient portfolios values:
  # Efficient frontier risk
  dbTable=data.frame(Date=tail(Returns$Date,1),
                     Value=portfolio1@portfolio@portfolio$targetRisk[,2],
                     Portfolio=experiment,
                     Ticker="Efficient frontier risk")
  # Efficient frontier return
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=portfolio1@portfolio@portfolio$targetReturn[,2],
                           Portfolio=experiment,
                           Ticker="Efficient frontier return") )
  # Efficient frontier CVaR
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=portfolio1@portfolio@portfolio$targetRisk[,3],
                           Portfolio=experiment,
                           Ticker="Efficient frontier CVaR") )
  # Efficient frontier VaR
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=portfolio1@portfolio@portfolio$targetRisk[,4],
                           Portfolio=experiment,
                           Ticker="Efficient frontier VaR") )
  # Efficient frontier covariance
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=portfolio1@portfolio@portfolio$targetRisk[,1],
                           Portfolio=experiment,
                           Ticker="Efficient frontier covariance") )
  
  # Efficient frontier weights:====
  
  weightsTable=cbind(data.frame(Date=tail(Returns$Date,1)),
                     portfolio1@portfolio@portfolio$weights*100,
                     Portfolio=seq(1:99))
  
  colnames(weightsTable)=c("Date",assetNames,"Portfolio")
  
  # Minimum variance portfolio estimation:
  print(paste0("Estimating minumum variance portfolio for portfolio: ",experiment))  
  mvportfolio=minvariancePortfolio(data=as.timeSeries(Returns),spec=portfoliospec)
  
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=mvportfolio@portfolio@portfolio$targetRisk[2],
                           Portfolio=experiment,
                           Ticker="Minimum variance portfolio risk") )
  
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=mvportfolio@portfolio@portfolio$targetReturn[2],
                           Portfolio=experiment,
                           Ticker="Minimum variance portfolio return") )
  
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=mvportfolio@portfolio@portfolio$targetRisk[3],
                           Portfolio=experiment,
                           Ticker="Minimum variance portfolio CVaR") )
  
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=mvportfolio@portfolio@portfolio$targetRisk[4],
                           Portfolio=experiment,
                           Ticker="Minimum variance portfolio VaR") )
  
  dbTable=rbind(dbTable,
                data.frame(Date=tail(Returns$Date,1),
                           Value=mvportfolio@portfolio@portfolio$targetRisk[1],
                           Portfolio=experiment,
                           Ticker="Minimum variance portfolio covariance") )
  # Min variance portfolio weights:
  weightsTable2=cbind(data.frame(Date=tail(Returns$Date,1)),
                      as.data.frame(t(mvportfolio@portfolio@portfolio$weights))*100,
                      data.frame(Portfolio="Minimum variance portfolio")
  )
  
  weightsTable=rbind(weightsTable,weightsTable2)
  
# outObject generation and call:====
  outObject=list(efficientFrontiers=dbTable,
                 portfolioWeights=weightsTable)
  return(outObject)
}

