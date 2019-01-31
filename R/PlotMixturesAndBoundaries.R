#' Shows GMM with Boundaries
#' 
#'  Plots Gaussian Mixture Model with Bayes decision boundaries, such that: \cr
#'  \cr
#'  Black is the PDE of Data \cr
#'  \cr
#'  Red is color of the GMM \cr
#'  \cr
#'  Magenta are the Bayes boundaries
#'
#' @param Data vector (1:N) of data points
#' @param Means vector[1:L] of Means of Gaussians (of GMM),L == Number of Gaussians
#' @param SDs  vector of standard deviations, estimated Gaussian Kernels, has to be the same length as Means
#' @param Weights vector of relative number of points in Gaussians (prior probabilities), has to be the same length as Means
#' @param IsLogDistribution Optional, ==1 if distribution(i) is a LogNormal, default vector of zeros of length 1:L
#' @param SingleColor Optional, Color for line plot of all the single gaussians, default magenta
#' @param MixtureColor Optional, Color of line plot for the mixture, default red
#' @param DataColor Optional, Color of line plot for the data, default black
#' @param BoundaryColor Optional, Color of bayesian boundaries
#' @param xlab Optional, x label, see plot
#' @param ylab Optional, y label, see plot
#' @param ... Optional, see plot for plot properties and for SingleGausses PlotMixtures
#'
#' @details 
#' you may also set SingleGausses=T than components of the mixture in blue will be shown.
#' @author Michael Thrun
#'
#'\strong{See Also}
#'
#' BayesDecisionBoundaries,PlotMixtures
#'

PlotMixturesAndBoundaries <-function(Data, Means, SDs, Weights, IsLogDistribution = rep(FALSE,length(Means)), SingleColor = 'blue', MixtureColor = 'red',DataColor='black', BoundaryColor = 'magenta', xlab, ylab, ...){
#PlotGaussMixesAndBoundaries(Data,Means,SDs,Weights,SingleColor,MixtureColor)
# Plot a Mixture of Gaussian/LogNormal and Bayesian decision boundaries

# INPUT
# Data[1:n]                   data column for which the distribution was modelled
# Means[1:L]                  Means of Gaussians,  L ==  Number of Gaussians
# SDs[1:L]                  estimated Gaussian Kernels = standard deviations
# Weights[1:L]                relative number of points in Gaussians (prior #															probabilities): sum(Weights) ==1
# OPTIONAL
# IsLogDistribution(1:L) 			gibt an ob die jeweilige Verteilung eine 
#															Lognormaverteilung ist,(default ==0*(1:L))
# SingleColor              		PlotSymbol of all the single gaussians, default 
#															magenta
# MixtureColor             		PlotSymbol of the mixture, default black
# BoundaryColor               PlotSymbol of the boundaries default green
# RoundNpower                 Decision Boundaries are rounded by 
#															roundn(DecisionBoundariesRoundNpower)
# ...													other plot arguments, like xlim = c(1,10)

# OUTPUT
# DecisonBoundaries(1:L-1)    Bayes decision boundaries 
# DBY(1:L-1)                  y values at the cross points of the Gaussians
       
# author MT 08/2015

# Labels
if(missing(xlab)){
	xlab = '' # no label for x axis
}
if(missing(ylab)){
	ylab = '' # no label for y axis
}

# Calculate intersections
dec = BayesDecisionBoundaries(Means,SDs,Weights,IsLogDistribution,Ycoor=T)

#MT: Bugifx
#if(is.list(dec)){
  DecisionBoundaries=as.vector(dec$DecisionBoundaries)
  #print('dec was a list, assuming usage of BayesDecisionBoundaries()')
#}
# Plot Gaussians
PlotMixtures(Data,Means,SDs,Weights,IsLogDistribution,SingleColor,MixtureColor,DataColor=DataColor, xlab = xlab, ylab  = ylab,...)

# intersecions
for (i in 1:length(DecisionBoundaries)){
 abline(v = DecisionBoundaries[i], col = BoundaryColor)
}

return (invisible(list(DecisionBoundaries = DecisionBoundaries, DBY = dec$DBY))) 
}

