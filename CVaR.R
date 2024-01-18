# CVaR function, given an estimated standard deviations vector and a confidence interval
# one

# Version 1.0 2024-01-18: Dr. Oscar Valdemar De la Torre Torres (in development for GED estimation).

library(fGarch)
library(doParallel)



CVaR=function(riskVector,confidenceVector,pdfFunct,CVaRt,Tlength){
library(fGarch)
cat("\f")
print("Setting CVaR parallel configurations up...")

  # parallel setup:
  #Setup backend to use many processors
  totalCores = detectCores()
  
  #Leave one core to avoid overload your computer
  cluster <- makeCluster(totalCores[1]-1) 
  registerDoParallel(cluster)
  
  cat("\f")
  print("Estimating CVaR...")
  
  CVaRTable=as.data.frame(matrix(0,length(riskVector),length(confidenceVector)))
  colnames(CVaRTable)=paste0("CVaR-",confidenceVector)
  
  for (a in 1:length(confidenceVector)){
   
 
    cat("\f")
    print(paste0("Estimating with ",pdfFunct,"pdf CVaR at ",confidenceVector[a]*100,"% of confidence..."))
    
# CVaR estimation===
    
    CVarVector=foreach(b=1:length(riskVector),.combine=c)%dopar%{  
#      parallelCVaR(riskVector[b],confidenceVector[a],pdfFunct,CVaRt,Tlength)
#---      
# parallel CVaR:
      alphaCVaR=1-confidenceVector[a]
      
      pValsSeq=seq(from=0,to=alphaCVaR,by=0.000001)
      pValsSeq=pValsSeq[-1]
      
      switch(pdfFunct,
             "norm"={
               zVal=qnorm(pValsSeq,0,1)
               cvar=mean((zVal*riskVector[b])*sqrt(CVaRt))
             },
             "t"={
               nu=Tlength-1
               tVal=qt(pValsSeq,nu)
               cvar=mean((tVal*riskVector[b])*sqrt(CVaRt))
             },
             "ged"={
               nu=1
               # q GED estimation:
               lambda = sqrt(2^(-2/nu) * gamma(1/nu)/gamma(3/nu))
               q = lambda * (2 * qgamma((abs(2 * pValsSeq - 1)), 1/nu))^(1/nu)
               gedVal = q * sign(2 * pValsSeq - 1) * 1 + 0
               
               cvar=mean((gedVal*riskVector[b])*sqrt(CVaRt))          
             }
      )
      
#---      
    } 

# CVaR table===      
    CVaRTable[,a]=CVarVector
     # a loop ends here: 
  }
  
  # output objects:
  return(CVaRTable)
  
  #Stopping parallel cluster
 stopCluster(cluster)  

  
}

