\name{dynPop}
\alias{estExp}
\alias{estDem}
\alias{crescLog}
\alias{popStr}
\alias{dynPop}

\title{
Population Dynamic Models
}

\description{
Functions to simulate population dynamic models. 
}

\usage{
estDem(tmax,n ,m ,N0 , nsim, ciclo)
estExp(N0,r,varr,tmax) 
crescLog(N0, r, K, tmax)
popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, rw, cl, tmax)
}

\arguments{
  \item{N0}{number of individuals at start time.
}
  \item{tmax}{maximum simulation time.
}
  \item{nsim}{number of populations simulated.
}
  \item{ciclo}{size of cicles in simulation.
}
  \item{cl}{number of columns to build the simulated scene.
}
  \item{rw}{number of rows to build the simulated scene.
}
  \item{r}{intrinsic growth rate.
}
  \item{varr}{variance.
}
  \item{n}{birth rate.
}
  \item{m}{dead rate.
}
  \item{K}{carrying capacity.
}
  \item{p.sj}{probability of seed survival.
}
  \item{p.jj}{probability of juvenil survival.
}
  \item{p.ja}{probability of transition from juvenil to adult fase.
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
}

\details{
'estExp' simulates an exponential population growth with environmental stochasticity. 

'estDem'simulates birth and death rates for stochastic (demographic) models.

Stochastic models uses 'r' and 'lambda' values taken of a normal distribution with mean 'r' and variance 'varr'.

'crescLog' simulates a logistic growth for continuous and discrete models.

'popStr' simulates a structured population dynamics, with Leftkovitch matrices.

In 'popStr' the number of patches in the simulated scene is defined by rw*cl.

}


\value{
The functions return graphics with the simulation results, and a matrix with the population size for deterministic and stochastic models.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\examples{
popStr(p.sj=0.4, p.jj=0.6, p.ja=0.2, p.aa=0.9, fec=0.8, ns=100,nj=40,na=20, rw=30, cl=30, tmax=20)
}

\keyword{population dynamics}
\keyword{simulation}
