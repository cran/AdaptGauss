#' BayesClassification
#' 
#' Bayes Klassifikation den Daten zuordnen
#' 
#'
#' @param Data vector of Data
#' @param Means vector[1:L] of Means of Gaussians (of GMM)
#' @param SDs vector of standard deviations, estimated Gaussian Kernels, has to be the same length as Means
#' @param Weights vector of relative number of points in Gaussians (prior probabilities), has to be the same length as Means
#' @param IsLogDistribution Optional, ==1 if distribution(i) is a LogNormal, default vector of zeros of length 1:L
#' @param ClassLabels Optional numbered class labels that are assigned to the classes. default (1:L), L number of different components of gaussian mixture model
#'
#' @return
#' Cls(1:n,1:d) classiffication of Data, such that 1= first component of gaussian mixture model, 2= second component of gaussian mixture model and so on. For Every datapoint a number is returned.
#' @author Michael Thrun
#' @export
#'
BayesClassification = function(Data,Means,SDs,Weights,IsLogDistribution=Means*0,
                               ClassLabels=c(1:length(Means))){
# Cls = BayesClassification(Data,M,S,W) 
# [Cls, DecisonBoundaries] = BayesClassification(Data,M,S,W,IsLogDistribution,ClassLabels) 
# Bayes Klassifikation den Daten zuordnen
# INPUT
# Data(1:n,1:d)             Data
# M(1:L),S(1:L),W(1:L)      parameters of the Gaussians/LogNormals
# 
# OPTIONAL
# IsLogDistribution(1:L)   ==1 if distribution(i) is a LogNormal, default Zeros
#
# OUTPUT
# Cls
# DecisonBoundaries(1:L-1)



DecisionBoundaries = BayesDecisionBoundaries(Means,SDs,Weights,IsLogDistribution) 
Cls = ClassifyByDecisionBoundaries(Data,DecisionBoundaries,ClassLabels) 

return(list(Cls=Cls, DecisionBoundaries=DecisionBoundaries))
}