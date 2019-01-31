#' Decision Boundaries calculated through Bayes Theorem
#'
#' Function finds the intersections of Gaussians or LogNormals
#'
#' @param Means vector[1:L] of Means of Gaussians (of GMM)
#' @param SDs vector of standard deviations, estimated Gaussian Kernels, has to be the same length as Means
#' @param Weights vector of relative number of points in Gaussians (prior probabilities), has to be the same length as Means
#' @param IsLogDistribution Optional, ==1 if distribution(i) is a LogNormal, default vector of zeros of length 1:L
#' @param MinData Optional, Beginning of range, where the Boundaries are searched for, default min(M
#' @param MaxData Optional, End of range, where the Boundaries are searched for, default max(M)
#' @param Ycoor Optional, Bool, if TRUE instead of vector of DecisionBoundaries list of DecisionBoundaries and DBY is returned
#'
#' @return
#' \describe{
#'   \item{DecisionBoundaries:}{vector[1:L-1], Bayes decision boundaries}
#'   \item{DBY:}{if (Ycoor==TRUE), y values at the cross points of the Gaussians is also returned, that the return is a list of DecisionBoundaries and DBY.}
#'   
#'}   
#' 
#' @author Michael Thrun, Rabea Griese
#'
#'@references
#'Duda, R. O., Hart, P. E., & Stork, D. G. (2001). Pattern classification. 2nd. Edition. New York, p. 512ff
#'
#'\strong{See Also}
#'
#' AdaptGauss,Intersect2Mixtures,Bayes4Mixtures
#' 
BayesDecisionBoundaries <- function(Means,SDs,Weights,IsLogDistribution,MinData,MaxData,Ycoor=F){
# DecisionBoundaries = BayesDecisionBoundaries(Means,SDs,Weights,IsLogDistribution);
# find the intersections of Gaussians or LogNormals
# 
# INPUT
# Means(1:L),SDs(1:L),Weights(1:L)      parameters of the Gaussians/LogNormals
# OPTIONAL
# IsLogDistribution(1:L)   ==1 if distribution(i) is a LogNormal, default zeros
# MinData
# MaxData
# Ycoor                   Bool, if TRUE instead of vector of DecisionBoundaries
#                         list of DecisionBoundaries and DBY is returned
# OUTPUT
# DecisionBoundaries(1:L-1)  Bayes decision boundaries 

# Author:  05/2015 RG
# 1.Editor: MT 08/2015
  
if(missing(IsLogDistribution)){
  IsLogDistribution = rep(FALSE,length(Means))
}  
if(missing(MinData)){
  MinData=(Means-3*SDs)[which.min(Means-3*SDs)]
}
if(missing(MaxData)){
  MaxData=(Means+3*SDs)[which.max(Means+3*SDs)]
}
L <- length(Means) # number of Gaussians

if(length(IsLogDistribution)!=L){
  warning(paste('Length of Means',L,'does not equal length of IsLogDistribution',length(IsLogDistribution),'Generating new IsLogDistribution'))  
  IsLogDistribution = rep(FALSE,L)
}
# sortieren nach Means
Sind = order(Means)
Means = sort(na.last=T,Means)
SDs = SDs[Sind]
Weights = Weights[Sind]
IsLogDistribution <- IsLogDistribution[Sind]

L1 = L-1

DecisionBoundaries <- Means[1:L1] # init
DBY = DecisionBoundaries*0 # init

# now intersect 2 Gauss
for(i in 1:L1){
	i1 = i+1
  Decision = Intersect2Mixtures(Means[i],SDs[i],Weights[i],Means[i1],SDs[i1],Weights[i1],IsLogDistribution[i:i1])
  DecisionBoundaries[i] = Decision$CutX
  DBY[i] = Decision$CutY
  DecisionBoundaries[i] = min(MaxData,DecisionBoundaries[i]);
  DecisionBoundaries[i] = max(MinData,DecisionBoundaries[i]);
}

if(Ycoor){
  return(list(DecisionBoundaries = DecisionBoundaries,DBY=DBY))
}else{
  return(DecisionBoundaries = DecisionBoundaries) 
}
}

