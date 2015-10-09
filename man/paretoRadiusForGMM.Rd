\name{paretoRadiusForGMM}
\alias{paretoRadiusForGMM}

\title{ParetoRadius for GMM}
\description{
  Calculation of the ParetoRadius i.e. the 18 percentiles of all mutual Euclidian distances in data.
}
\usage{
paretoRadiusForGMM(Data,maximumNrSamples,plotDistancePercentiles)
}
\arguments{
  \item{Data}{numeric data vector}
  \item{maximumNrSamples}{Optional, numeric. Maximum number for which the distance calculation can be done. 1000 by default.}
  \item{plotDistancePercentiles}{Optional, logical. If TRUE,  a plot of the percentiles of distances is produced. FALSE by default.}
}

\value{
  numeric value, the Pareto radius.
}
\details{The Pareto-radius of the hyperspheres is derived from the optimation of information for minimal set size.}
\note{
  paretoRadiusForGMM() is a kernel density estimation for variable intervals.
}
\references{Ultsch, A.: Pareto Density Estimation: A Density Estimation for Knowledge Discovery, Baier D., Wernecke K.D. (Eds), In Innovations in Classification, Data Science, and Information Systems - Proceedings 27th Annual Conference of the German Classification Society (GfKL) 2003, Berlin, Heidelberg, Springer, pp, 91-100, 2003.}
\author{Michael Thrun}
\seealso{\code{\link{paretoDensityEstimationForGMM}}, \code{\link{OptimalNoBins}}}
\keyword{kernel density estimation}