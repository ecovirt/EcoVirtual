\name{dynPop}
\alias{estExp}
\alias{estDem}
\alias{crescLog}
\alias{popStr}

\title{
Population Dynamic Models
}

\description{
Functions to simulate population dynamic models. 
}

\usage{
estDem(N0, b, d, tmax)
estExp(N0,r,varr,tmax) 
crescLog(N0, r, K, tmax)
popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, rw, cl, tmax)
}

\arguments{
  \item{N0}{number of individuos at start time.
}
  \item{tmax}{maximum time.
}
  \item{cl}{number of columns.
}
  \item{rw}{number of rowns.
}
  \item{r}{intrinsic growth rate.
}
  \item{varr}{variance.
}
  \item{b}{birth rate.
}
  \item{d}{dead rate.
}
  \item{K}{carrying capacity.
}
  \item{p.sj}{probability of seed survival.
}
  \item{p.jj}{probability of juvenil survival.
}
  \item{p.ja}{juvenil probability of transition to adult fase.
}
  \item{p.aa}{probability of adult survival.
}
  \item{ns}{number of seeds at initial time.
}
  \item{na}{number of adults at initial time.
}
  \item{nj}{number of juvenil at initial time.
}
  \item{fec}{mean number of propagulus per adult each cicle.
}

\details{
'estExp' gives the results for exponential population growth with environmental stochasticity in time 

'estDem' gives the results for the birth and death rates for deterministic and stochastic (enviromental) models in time.

Stochastic models uses 'r' and 'lambda' values taken of a normal distribution with mean 'r' and variance 'varr'.

'crescLog' gives the results for a logistic growth for continuous and discrete models in time.

'popStr' gives the result for a structured population dynamics, with Leslie's matrices for the stages of the population.


 
}


\value{
The functions return graphics with the simulation results, and an invisible matrix with the population size for deterministc and stochastic models on each time.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\examples{
popStr(p.sj=0.4, p.jj=0.6, p.ja=0.2, p.aa=0.9, fec=0.8, ns=100,nj=40,na=20, ln=30, cl=30, tmax=20)
}

\keyword{population dynamics}
\keyword{simulation}
