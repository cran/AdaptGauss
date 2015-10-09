\name{Chi2testMixtures}
\alias{Chi2testMixtures}

\title{Pearson's chi-squared test}
\description{
  Returns a P value and visualizes for chi-square test of Data versus a given Gauss Mixture Model
}
\usage{
Chi2testMixtures(Data,Means,SDs,Weights,IsLogDistribution,PlotIt,UpperLimit,VarName)
}

\arguments{
  \item{Data}{vector of data points}
  \item{Means}{vector of Means of Gaussians}
  \item{SDs}{vector of standard deviations, estimated Gaussian Kernels}
  \item{Weights}{vector of relative number of points in Gaussians (prior probabilities)}
  \item{IsLogDistribution}{Optional, if IsLogDistribution(i)==1, then mixture is lognormal, default vector of zeros of length 1:L}
  \item{PlotIt}{Optional, Default: FALSE, do a Plot of the compared cdfs and the KS-test distribution (Diff)}
  \item{UpperLimit}{Optional. test only for Data <= UpperLimit, Default = max(Data) i.e all Data.}
	\item{VarName}{If PlotIt=TRUE, the name of the inspected variable, default 'Data'}
}
\details{
  Let O_i be the observed features and E_i be the expected number E,
  than the test statistic is defined with T=sum((O_i-E_i)^2/E_i)*1/m, where m the number of data points.
  The expected number Ei may be derived for each bin with the null hypothesis that the data is from the distribution.
}
\value{
List with
  \item{Pvalue }{Pvalue of a suiting chi-square , Pvalue ==0 if Pvalue <0.001}
  \item{BinCenters }{bin centers}
  \item{ObsNrInBin }{No. of data in bin}
	\item{ExpectedNrInBin }{No. of data that should be in bin according to GMM}
	\item{Chi2Value }{the  TestStatistic i.e.: sum((ObsNrInBin(Ind)-ExpectedNrInBin(Ind))^2/ExpectedNrInBin(Ind)) with
                                 Ind = find(ExpectedNrInBin>=10)}
}
\references{Hartung, J., Elpelt, B., and Kloesener, K.H.: Statistik, 8. Aufl. Verlag Oldenburg (1991).}
\author{Rabea Griese}

