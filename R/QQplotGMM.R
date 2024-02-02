QQplotGMM=function(Data, Means, SDs, Weights, IsLogDistribution = Means*0,
                   Method = 1,
                   Line = TRUE, PlotSymbol = 20, col = "red",
                   xug = NULL, xog = NULL, LineWidth = 2, PointWidth = 0.8,
                   PositiveData = FALSE, Type = 8, NoQuantiles = 10000,
                   ylabel = 'Data', main = 'QQ-plot Data vs GMM',
                   lwd = 3, pch = 20, xlabel='Gaussian Mixture Model', ...){
  # QQplotGMM(Data,Means,SDs,Weights,IsLogDistribution,Line,PlotSymbol,xug,xog,LineWidth,PointWidth)
  # Quantile/Quantile = QQ-Plot im Vergleich. zu einem Gauss Mixture Model oder LGL Model
  # INPUT
  # Data(1:n)	             Daten, deren Verteilung verglichen werden soll
  # Means(1:L), SDs(1:L), Weights(1:L) die Paramter von Gaussians N(i) = Weights(i) * N(Means(i),SDs(i)
  #                        die Gesamtverteilung ergibst sich als Summe der N(i)
  # OPTIONAL
  # Method                 Integer: Old version stays as default (Method = 1)
  #                                 New method: Method == 2 (enforces new properties and robustness)
  # IsLogDistribution(1:L) gibt an ob die Einzelverteilung einer (generalisierten)Lognormaverteilung ist
  #                        wenn IsLogDistribution(i)==0 dann Mix(i) = Weights(i) * N(Means(i),SDs(i)
  #                        wenn IsLogDistribution(i)==1 dann Mix(i) = Weights(i) * LogNormal(Means(i),SDs(i)
  #                        Default: IsLogDistribution = Means*0;
  # Line									Line in QQplot: =TRUE (Default), without False
  # PlotSymbol             Symbol fur den qqplot, wenn nicht angegeben: PlotSymbol='b.'
  # col                    Character: color of regression line (only for Method = 2)
  # Type                   Integer: number of method used for computing the quantiles
  # NoQuantiles            Integer: Number of quantiles to compute (only for Method = 2)
  # xug,xog                Grenzen der Interpolationsgeraden,  interpoliert wird fuer percentiles(x) in [xug,xog]
  #                        Default: xug==min(x),xog==max(x), MT: Noch nicht implementiert!
  # LineWidth              Linienbreite der Interpolationsgeraden; Default =3
  # PointWidth             Dicke der Punkte im QQplot, existert nicht in Matlab
  # LineSymbol             Liniensymbol  der Interpolationsgerade;  Default ='r-'   MT: Nicht Implementiert
  # lwd                    Integer: graphic parameter - line width option (only for Method = 2)
  # pch                    Integer: graphic parameter for points (only for Method = 2)
  # in \dbt\Plot
  
  # benutzt randomLogMix und qqplotfit
  # MT 2014, reimplementiert aus Matlab von ALU 
  # Aus historischen Gr?nden QQplotGMM MIT Ausgleichgerade
  # QMS 2023: Integrate the DataVisualization approach of QQPlot (Method == 2)
  
  #xug = min(Data);  xog = max(Data); zu implementieren
  # LineSymbol='r-' nicht implementiert
  
  
  GMM = RandomLogGMM(Means,SDs,Weights,IsLogDistribution);
  if(PositiveData == TRUE){
    GMM = GMM[GMM >= 0]
  }
  
  if(Method == 1){
    quants<-qqplot(GMM, Data, pch=PlotSymbol, col="blue", cex=PointWidth, xlab=xlabel, ylab=ylabel, main=main,...) #MT: na.rm=TRUE argument weglassen
    if(Line){
      fit<-lm(quants$y~quants$x)
      summary(fit)
      abline(fit, col="red", lwd=LineWidth)
    }
    return(invisible(quants))
  }else{
    Res = DataVisualizations::QQplot(X = GMM, Y = Data, Type = Type,
                               NoQuantiles = NoQuantiles,
                               xlab = xlabel, ylab = ylabel, col = col,
                               main = main, lwd = lwd, pch = pch, ...)
    return(Res)
  }
}