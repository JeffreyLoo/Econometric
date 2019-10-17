library(fBasics)
library(rugarch) 
library(forecast)
library(sos)
library(quantmod)

##################### P A R T A #######################

###Opening data files ###
btc = read.table("btc.txt", header=T) # Load the data
eth = read.table("eth.txt", header=T) # Load the data


###Calculating the log difference ###
logbtc =as.ts(100*diff(log(btc$btc)))
logeth =as.ts(100*diff(log(eth$eth)))



##################### P A R T B #######################

###This function fits and predicts the one day ahead forecast for each rolling sample for the 8 models
###parameter data: the original log returns of cyptocurrency
###parameter size: the size of the rolling window sample
###loop A is in charge of making 8 model for each rolling windows e.g. making 8 models in 1 shot for one rolling window sample
###Models are, in sequence:
### 1.AR(1)-GARCH(1,1) 
### 2.AR(1)-GJR-GARCH(1,1) 
### 3.AR(5) 
### 4.Naive 
### 5.Historical Mean
### 6.Simple Moving Average(20)
### 7.Simple Moving Average(60)
### 8.Simple Moving Average(120)

predicted <- function(data, size ){
  
  #holds the predictions for 8 models
  #each models will have 500 observations
  #so result will have 500*8 observations
  result =list(c(),c(),c(),c(),c(),c(),c(),c())
  
  #to get r^2
  datasq =sapply(data, function(x) x^2)
  
  #the specification for sGARCH model
  spec1 <- ugarchspec(variance.model = list(model = "sGARCH", 
                                            garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), arfima = FALSE, 
                                        include.mean = TRUE),
                      distribution.model = "norm") 
  
  #the specification for gjrGARCH model
  spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), arfima = FALSE, include.mean = TRUE),
                      distribution.model = "norm")
  
  #A
  for (i in seq(1,500,1)){ 
    
    #AR(1)-GARCH(1,1)
    sgarch.fit = ugarchfit(data[i:(size+i)], spec = spec1, solver = "hybrid") # Model estimation 
    forc1 = ugarchforecast(sgarch.fit, n.ahead=1)
    result[[1]][[i]] =attributes(forc1)[[1]]$sigmaFor[1]^2
    
    
    #AR(1)-GJR-GARCH(1,1) 
    gjrgarch.fit =ugarchfit( data[i:(size+i)], spec = spec2, solver = "hybrid")
    forc2 = ugarchforecast(gjrgarch.fit, n.ahead=1)
    result[[2]][[i]] =attributes(forc2)[[1]]$sigmaFor[1]^2
    
    #AR(5) model
    ARmodel =arima(datasq[(i):(size+i)], order=c(5,0,0), method="ML", optim.method="BFGS")
    result[[3]][[i]] =predict(ARmodel,1)$pred[1]
    
    #Naive model
    result[[4]][[i]] =datasq[size+i]
    
    #Historical Mean Approach
    result[[5]][[i]] =mean(datasq[(i):(size+i)])
    
    #Simple Moving Average
    counter = 6
    
    #B this loop is to calculate the simple moving average 20,60,180
    for (j in c(20,60,180)){
      result[[counter]][[i]] =mean(datasq[(size+i-j+1):(size+i)])
      counter =counter +1
      
    }
  }
  return (result)
}







### This function calculates the MSE, MAD, and QLIKE 
###parameter predict: the predicted values from 8 models
###parameter data: the last 500 observations of the original log return of cryptocurrency
### The results are separated in 3 different loss functions.
### each loss functions will have ordered calculation from 8 models to determine which model has the lowest loss function value

accurateModel <- function(predict,data){
  
  #to square the data
  datasq =sapply(data, function(x) x^2)
  
  #the table to store transformed observations of the predictions
  result=c()
  
  #store 3 tables each table holds value of loss functions in order
  error=c()
  
 
  #MSE
  result[[1]] =cbind(cbind((predict[[1]]-datasq)^2),
                     cbind((predict[[2]]-datasq)^2),
                     cbind((predict[[3]]-datasq)^2),
                     cbind((predict[[4]]-datasq)^2),
                     cbind((predict[[5]]-datasq)^2),
                     cbind((predict[[6]]-datasq)^2),
                     cbind((predict[[7]]-datasq)^2),
                     cbind((predict[[8]]-datasq)^2))
  
  #MAD
  result[[2]] =cbind(cbind(abs(predict[[1]]-datasq)),
                     cbind(abs(predict[[2]]-datasq)),
                     cbind(abs(predict[[3]]-datasq)),
                     cbind(abs(predict[[4]]-datasq)),
                     cbind(abs(predict[[5]]-datasq))
                     ,cbind(abs(predict[[6]]-datasq)),
                     cbind(abs(predict[[7]]-datasq)),
                     cbind(abs(predict[[8]]-datasq)))
  
  #QLIKE
  result[[3]] =cbind(cbind(log((predict[[1]])^2)+((predict[[1]])^-2)*datasq^2),
                     cbind(log((predict[[2]])^2)+((predict[[2]])^-2)*datasq^2),
                     cbind(log((predict[[3]])^2)+((predict[[3]])^-2)*datasq^2),
                     cbind(log((predict[[4]])^2)+((predict[[4]])^-2)*datasq^2),
                     cbind(log((predict[[5]])^2)+((predict[[5]])^-2)*datasq^2),
                     cbind(log((predict[[6]])^2)+((predict[[6]])^-2)*datasq^2),
                     cbind(log((predict[[7]])^2)+((predict[[7]])^-2)*datasq^2),
                     cbind(log((predict[[8]])^2)+((predict[[8]])^-2)*datasq^2))

  #to name the columns
  #to find average of transformed prediction observations to obtain a single value of the loss function
  for (i in 1:3){
    colnames(result[[i]]) <- c("ARGARCH(1,1)","ARGJRGARCH(1,1)","AR(5)","Naive","Historical Mean","SMA(20)","SMA(60)","SMA(180)")
    mean <- colMeans(result[[i]])
    error[[i]] <- cbind(mean)
    error[[i]] <- error[[i]][order(as.numeric(error[[i]][,1])), ]
    
  }
  

  
  return (error)
  
}





# (i)
###calling functions here

###Predicted results for 8 models
predicted_btc=predicted(logbtc,1502) 
predicted_eth=predicted(logeth,976) 



###Error for each models for each cryptocurrencies in ascending order
error_btc =accurateModel(predicted_btc,tail(logbtc,500)) 
error_eth =accurateModel(predicted_eth,tail(logeth,500)) 



#set the display to print more
options("max.print" =100000)

