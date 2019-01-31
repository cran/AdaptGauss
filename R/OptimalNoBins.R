#' Optimal Number Of Bins
#' 
#' Calculation of the optimal number of bins for a histogram.
#'
#' @param Data Data
#' 
#' @details 
#' 
#' The bin width ist defined with bw=3.49*stdrobust(1/(n)^1/3)
#'
#' @return
#' 
#' optNrOfBins The best possible number of bins. Not less than 10 though
#' 
#' @note 
#' OptimalNoBins() is a kernel density estimation for fixed intervals.
#' 
#' @author Alfred Ultsch, Michael Thrun
#' 
#' @references 
#' 
#' David W. Scott Jerome P. Keating: A Primer on Density Estimation for the Great Home Run Race of 98, STATS 25, 1999, pp 16-22.
#'
#' \strong{See Also}
#' 
#' ParetoRadius
#' 
#' @examples
#' 
#' Data = c(rnorm(1000),rnorm(2000)+2,rnorm(1000)*2-1)
#' 
#' optNrOfBins = OptimalNoBins(Data)
#' 
#' minData = min(Data,na.rm = TRUE)
#' 
#' maxData = max(Data,na.rm = TRUE)
#' 
#' i = maxData-minData
#' 
#' optBreaks = seq(minData, maxData, i/optNrOfBins) # bins in fixed intervals 
#' 
#' hist(Data, breaks=optBreaks)
#' 
OptimalNoBins <-
function(Data){
# function [OptimalNrOfBins] = OptNrOfBins(data); 
# % [OptimalNrOfBins] = OptNrOfBins(data)
# %
# % DESCRIPTION
# % Berechung der optimalen Anzahl von Bins fuer ein Histogramm
# % nach Keating/Scott 99
# % INPUT
# % data               die Daten
# % OUTPUT
# % OptimalNrOfBins   die bestmoegliche ANzahl von Bins, minimal jedoch 10
# %                   Verwendung fuer hist(data,OptimalNrOfBins);
  #Anzahl vorhandene Daten
    if(is.matrix(Data)) nData <- colSums(!is.nan(Data))
    if(is.vector(Data)) nData <- sum(!is.nan(Data))
    
    prctile_hlp=function(x,p){
#   matlab:
#   Y = prctile(X,p) returns percentiles of the values in X. 
#   p is a scalar or a vector of percent values. When X is a 
#   vector, Y is the same size as p and Y(i) contains the p(i)th 
#   percentile. When X is a matrix, the ith row of Y contains the 
#   p(i)th percentiles of each column of X. For N-dimensional arrays,
#   prctile operates along the first nonsingleton dimension of X.  
  if(length(p)==1){  
            if(p>1){p=p/100}
            
  }

  if(is.matrix(x) && ncol(x)>1){
    cols<-ncol(x)
    quants<-matrix(0,nrow=length(p),ncol=cols)
    for(i in 1:cols){
      quants[,i]<-quantile(x[,i],probs=p,type=5,na.rm=TRUE)
    }
  }else{
    quants<-quantile(x,p,type=5,na.rm=TRUE)
  }
  return(quants)
} 
    
    
    if(nData<1){
      optNrOfBins<-0
    }else{    
      sigma<-sd(Data,na.rm=TRUE)    
      p<-prctile_hlp(Data,c(0.25,0.75))
      interquartilRange<-p[2]-p[1]     
      sigmaSir<-min(sigma,interquartilRange/1.349)
      optBinWidth<-3.49*sigmaSir/(nData)^(1/3)
      if(optBinWidth>0){
        optNrOfBins<-max(ceiling((max(Data,na.rm=TRUE)-min(Data,na.rm=TRUE))/optBinWidth),10)
      }else{
        optNrOfBins<-10
      }
    }                  
 return (optNrOfBins) 
 }