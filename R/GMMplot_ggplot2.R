GMMplot_ggplot2 <- function(Data, Means, SDs, Weights, BayesBoundaries, SingleGausses = TRUE, Hist = FALSE)
{
  
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    message(
      'ggplot2 package is needed for plotting. Please install the package which is defined in "Suggests".') 
    return(NULL)
  }
  
  if (!requireNamespace('reshape2', quietly = TRUE)) {
    message(
      'reshape2 package is missing. Please install the package which is defined in "Suggests".') 
    return(NULL)
  }
  
  
  if(missing(BayesBoundaries))
    BayesBoundaries=AdaptGauss::BayesDecisionBoundaries(Means, SDs, Weights,Ycoor = FALSE)
  
  MeansO <- Means[order(Means)]
  SDsO <- SDs[order(Means)]
  WeightsO <- Weights[order(Means)]
  PDE <- DataVisualizations::ParetoDensityEstimation(Data)
  
  GMMn <- data.frame(cbind(PDE$kernels, PDE$paretoDensity))
  for (i in 1:length(MeansO)) GMMn <- cbind(GMMn, WeightsO[i]*dnorm(PDE$kernels, mean=MeansO[i],sd=SDsO[i],log=F))
  GMMn <- cbind(GMMn, rowSums(GMMn[3:(2+length(MeansO))]))
  names(GMMn) <- c("kernels", "density", paste0("M",1:i), "Sum")
  dfGMM <- reshape2::melt(data = GMMn, id="kernels")
  dfGMM$Cls=factor(dfGMM$variable)
  
  #require(grDevices)
  if (SingleGausses == TRUE & length(Means) > 1) p1 <- ggplot2::ggplot(data = dfGMM, ggplot2::aes_string(colour = 'Cls')) else 
    p1 <- ggplot2::ggplot(data = subset(dfGMM, dfGMM$variable %in% c("density", "Sum")), ggplot2::aes_string(colour = 'Cls'))
  if(Hist == TRUE) {
    breaks <- pretty(range(Data), n = nclass.FD(Data), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    p1 <- p1 + ggplot2::geom_histogram(data = data.frame(Data), ggplot2::aes_string(x = Data, "..density.."), fill = "grey85", color = "grey95", binwidth=bwidth) + ggplot2::guides(fill = FALSE) 
  }
  GMMplot <- p1 + ggplot2::geom_line(ggplot2::aes_string(x = 'kernels', y = 'value'),  size = ggplot2::rel(0.9)) + ggplot2::labs(colour = "Lines") + ggplot2::labs(x = "Data", y = "Probability Density")
  if(length(MeansO) > 1) GMMplot <- GMMplot + ggplot2::geom_vline(xintercept=BayesBoundaries, color = "magenta", linetype = 2) +
    ggplot2::annotate("text", x = MeansO, y = rep(.01,length(MeansO)), size = ggplot2::rel(4), label = paste0("M#", 1:length(MeansO)))
  return(GMMplot)
}
