\name{LikelihoodRatio4Mixtures}
\alias{LikelihoodRatio4Mixtures}

\title{Likelihood Ratio for Gaussian Mixtures}
\description{
  Computes the likelihood ratio for two Gaussian Mixture Models.
}
\usage{
LikelihoodRatio4Mixtures(Data,NullMixture,OneMixture,PlotIt,LowerLimit,UpperLimit)
}

\arguments{
  \item{Data}{ Data points.}
  \item{NullMixture}{ A Matrix: cbind(Means0,SDs0,Weights0) or  cbind(Means0,SDs0,Weights0,IsLog0). 
	The null model; usually with less Gaussians than the OneMixture}
  \item{OneMixture}{ A Matrix: cbind(Means1,SDs1,Weights1) or  cbind(Means1,SDs1,Weights1,IsLog1). 
	The alternative model usually with more Gaussians than the OneMixture.}
  \item{PlotIt}{ Optional: zero or one. o a Plot of the compared cdf's and the KS-test distribution (Diff)}
  \item{LowerLimit}{ Optional: test only for Data >= LowerLimit, Default = min(Data) i.e all Data.}
  \item{UpperLimit}{ Optional: test only for Data <= UpperLimit, Default = max(Data) i.e all Data.}
  
	}
\value{
List with
  \item{Pvalue}{the error that we make, if we accept OneMixture as the better Model over the NullMixture}
  \item{NullLogLikelihood}{log likelihood of GMM Null}
  \item{OneLogLikelihood}{log likelihood of GMM One}
}

\author{Alfred Ultsch, Michael Thrun, Catharina Lippmann}

