\name{DistributionOptimization}
\alias{DistributionOptimization}

\title{DistributionOptimization}
\description{This function models a mixture of gaussian distributions by minimizing the chi2 value for a data sample. An evolutional algorithm is used to achieve that.}

\arguments{
\item{Data[1:n]}{Data to be modelled}
\item{Modes}{Number of expected Modes}
\item{Monitor=1}{0: no monitoring, 1: status messages, 2: status messages and plots}
\item{SelectionMethod=4}{ 1: linear rank selection
     2: non linear rank selection
     3: proportional roulette wheel
     4: unbiased tournament
     5: fitness proportional selection
     6: fitness proportional mit goldberg sigma truncation
     7: biased tournament}
\item{MutationMethod=6}{ 1: uniform random mutation
     2: nonuniform random mutation
     3: random mutation around solution
     4: GaussMutationWithIntensity,  alternative random mutation around solution
     5: GaussMutationCust
     6: TwoPhaseMutation - mutation is uniform random during the first half of iterations, and than focuses around current solution}
\item{CrossoverMethod=2}{1: single point crossover
     2: whole arithmetic crossover
     3: local arithmetic crossover
     4: blend crossover
     5: GaussCrossover - exchange complete gausscomponents
     6: MultiPointCrossover  - Random amount of information between mixtures get exchanged}
\item{PopulationSize=Modes*60}{size of the population}
\item{Mutationrate=0.5}{amount (0..1) of population that gets mutated}
\item{SortedGauss=T}{should the gauss components be sorted by expactation value after each generation?}
\item{Elitism=0.05}{amount of best individuals that will survive generation unchanged}
\item{CrossoverRate=0.5}{amount of individuals that will be used for crossover}
\item{Iter=Modes*500}{number of iterations of this algorithm}
\item{OverlapRatio=0.8}{ratio between Chi2 and OverlapError (only if FitnessMethod = 4 (Chi2ValueWithOverlap))}
\item{IsLogDistribution=rep(F, Modes)}{which gauss components should be considered as log gaussian}
\item{Means=NULL}{startpopulation is initialized with this Means}
\item{SDs=NULL}{startpopulation is initialized with this SDs}
\item{Weights=NULL}{startpopulation is initialized with this Weights}

}
\value{
  List with
  \item{GA}{GA object of the training}
  \item{Weights}{...}
  \item{SDs}{...}
  \item{Means}{...}
  \item{precalc}{precalc object (for debugging purposes)}
}
