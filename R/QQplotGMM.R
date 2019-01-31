#' Quantile Quantile Plot of Data
#' 
#' Quantile Quantile plot of data against gaussian distribution mixture model with optional best-fit-line
#'
#'
#' @param Data vector (1:N) of data points
#' @param Means vector[1:L] of Means of Gaussians (of GMM),L == Number of Gaussians
#' @param SDs vector of standard deviations, estimated Gaussian Kernels, has to be the same length as Means
#' @param Weights vector of relative number of points in Gaussians (prior probabilities), has to be the same length as Means
#' @param IsLogDistribution Optional, ==1 if distribution(i) is a LogNormal, default Zeros of Length L
#' @param Line Optional, Default: TRUE=Regression Line is drawn
#' @param PlotSymbol Optional, Symbol to be used within qqplot
#' @param xug Optional, lower limit of the interval [xug, xog], in which a line will be interpolated
#' @param xog Optional, upper limit of the interval [xug, xog], in which a line will be interpolated
#' @param LineWidth Optional, width of regression line, if Line==TRUE
#' @param PointWidth Optional, width of points
#' @param ylab Optional, see plot
#' @param main Optional, see plot
#' @param ... Note: xlab cannot be changed, other parameters see qqplot
#' 
#' 
#' @details 
#' Only verified for a Gaussian Mixture Model, usage of IsLogDistribution for LogNormal Modes is experimental!
#'
#' @return List With 
#' \describe{
#'   \item{x:}{The x coordinates of the points that were plotted}
#'   \item{y:}{The original data vector, i.e., the corresponding y coordinates}
#'}
#'
#' @author Michael Thrun 
#' 
#' @references 
#' Michael, J. R. (1983). The stabilized probability plot. Biometrika, 70(1), 11-17.
#' 
#' \strong{See Also}
#' 
#' qqplot
#'
#' @examples
#' 
#' data=c(rnorm(1000),rnorm(2000)+2,rnorm(1000)*2-1)
#' QQplotGMM(data,c(-1,0,2),c(2,1,1),c(0.25,0.25,0.5))
QQplotGMM=function(Data,Means,SDs,Weights,IsLogDistribution=Means*0,Line=TRUE,PlotSymbol=20,
                   xug=NULL,xog=NULL,LineWidth=2,PointWidth=0.8,
                   ylab='Data',main='QQ-plot Data vs GMM',...)
# QQplotGMM(Data,Means,SDs,Weights,IsLogDistribution,Line,PlotSymbol,xug,xog,LineWidth,PointWidth)
# Quantile/Quantile = QQ-Plot im Vergleich. zu einem Gauss Mixture Model oder LGL Model
# INPUT
# Data(1:n)	             Daten, deren Verteilung verglichen werden soll
# Means(1:L), SDs(1:L), Weights(1:L) die Paramter von Gaussians N(i) = Weights(i) * N(Means(i),SDs(i)
#                        die Gesamtverteilung ergibst sich als Summe der N(i)
# OPTIONAL
# IsLogDistribution(1:L) gibt an ob die Einzelverteilung einer (generalisierten)Lognormaverteilung ist
#                        wenn IsLogDistribution(i)==0 dann Mix(i) = Weights(i) * N(Means(i),SDs(i)
#                        wenn IsLogDistribution(i)==1 dann Mix(i) = Weights(i) * LogNormal(Means(i),SDs(i)
#                        Default: IsLogDistribution = Means*0;
# Line									Line in QQplot: =TRUE (Default), without False
# PlotSymbol             Symbol fur den qqplot, wenn nicht angegeben: PlotSymbol='b.'
# xug,xog                Grenzen der Interpolationsgeraden,  interpoliert wird fuer percentiles(x) in [xug,xog]
#                        Default: xug==min(x),xog==max(x), MT: Noch nicht implementiert!
# LineWidth              Linienbreite der Interpolationsgeraden; Default =3
# PointWidth             Dicke der Punkte im QQplot, existert nicht in Matlab
# LineSymbol             Liniensybol  der Interpolationsgerade;  Default ='r-'   MT: Nicht Implementiert
#
# in \dbt\Plot

# benutzt randomLogMix und qqplotfit
# MT 2014, reimplementiert aus Matlab von ALU 
# Aus historischen Gr?nden QQplotGMM MIT Ausgleichgerade

{
  

#xug = min(Data);  xog = max(Data); zu implementieren
# LineSymbol='r-' nicht implementiert
xlabel='Gaussian Mixture Model'

GMM = RandomLogGMM(Means,SDs,Weights,IsLogDistribution);

 quants<-qqplot(GMM, Data, pch=PlotSymbol, col="blue", cex=PointWidth, xlab=xlabel, ylab=ylab, main=main,...) #MT: na.rm=TRUE argument weglassen
 if(Line){
 fit<-lm(quants$y~quants$x)
 summary(fit)
 abline(fit, col="red", lwd=LineWidth)
 }
 return(invisible(quants))
}