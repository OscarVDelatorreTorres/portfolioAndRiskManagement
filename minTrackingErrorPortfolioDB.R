# This is the min tracking error portfolio function for the ICTI 2023-2024 research proyect of
# Oscar De la Torre-Torres <https://www.oscardelatorretorres.com/icti2023>
# This routine estimates the minimum tracking error portfolio for a given benchmark portfolion
# In this case non-commodities such as avocado or green lemon.
# The idea of this function is to estimate the optimal portfolio results, along with the
# parameters of the minimum variance portfolio (minimum tracking error).
# The core idea is to estimate, given a Returns data matrix and a preselected
# portfolio parameters function, the theoretical portfolios and store the results in a database.

minTrackingErrorPortfolio=function(Returns,myestimator,connDB,experiment){
  
  outObject=list()
  return(outObject)
}


#Portfolio frontiers======  
# Conventional portfolio spec:
portfoliospec=portfolioSpec(portfolio= list(
  weights = NULL, targetReturn = NULL,
  targetRisk = NULL, riskFreeRate = 0, nFrontierPoints = 100,
  status = NA))

setEstimator(portfoliospec)=myestimator

# Efficient frontiers estimation
portfolio1=portfolioFrontier(data=as.timeSeries(Returns),spec=portfoliospec)

dbTable=data.frame(Risk=portfolio1@portfolio@portfolio$targetRisk[,2],
                     Return=portfolio1@portfolio@portfolio$targetReturn[,2],
                     Portfolio=experiment,
                   Ticker="Efficient frontier")

mvportfolio=minvariancePortfolio(data=as.timeSeries(Returns),spec=portfoliospec)

