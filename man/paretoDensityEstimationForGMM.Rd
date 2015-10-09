\name{paretoDensityEstimationForGMM}
\alias{paretoDensityEstimationForGMM}
\title{Pareto Density Estimation}
\description{
This function estimates the Pareto Density for the distribution of one variable.}
\usage{
paretoDensityEstimationForGMM(Data,paretoRadius,kernels,MinAnzKernels)
}
\arguments{
  \item{Data}{numeric vector of data.}
  \item{paretoRadius}{Optional, numeric value, see  \code{\link{paretoRadiusForGMM}}}
  \item{kernels}{Optional, numeric vector. data values where pareto density is measured at. If 0 (by default) kernels will be computed. }
    \item{MinAnzKernels}{Optional, minimal number of kernels, default MinAnzKernels==100 }
}
\value{

  \item{kernels }{numeric vector. data values at with Pareto Density is measured.}
  \item{paretoDensity}{numeric vector containing the determined density by ParetoRadius.}
  \item{paretoRadius}{numeric value.}
}
\details{
Pareto Density Estimation (PDE) is a method for the estimation of probability density functions using hyperspheres. The Pareto-radius of the hyperspheres is derived from the optimation of information for minimal set size. It is shown, that Pareto Density is the best estimate for clusters of Gaussian structure. The method is shown to be robust when cluster overlap and when the variances differ across clusters.
}
\author{Michael Thrun}
\note{This is the best density estimation to judge Gaussian Mixtures  of the data see [Ultsch 2003]
}
\examples{
			data = c(rnorm(1000),rnorm(2000)+2,rnorm(1000)*2-1)
			# Optional, if not given paretoDensityEstimationForGMM() 
			# calls paretoRadiusForGMM()
			pareto_radius<-paretoRadiusForGMM(data) 
			pdeVal        <- paretoDensityEstimationForGMM(data,pareto_radius)
			plot(pdeVal$kernels,pdeVal$paretoDensity,type='l',xaxs='i',
			yaxs='i',xlab='Data',ylab='PDE') 
  			
}
\references{Ultsch, A.: Pareto Density Estimation: A Density Estimation for Knowledge Discovery, Baier D., Wernecke K.D. (Eds), In Innovations in Classification, Data Science, and Information Systems - Proceedings 27th Annual Conference of the German Classification Society (GfKL) 2003, Berlin, Heidelberg, Springer, pp, 91-100, 2003.}
\seealso{\code{\link{paretoRadiusForGMM}}}