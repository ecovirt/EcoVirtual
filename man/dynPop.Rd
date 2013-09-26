\name{dynPop}
\alias{popExp}
\alias{estEnv}
\alias{estDem}
\alias{popLog}
\alias{popStr}
\alias{dynPop}
\alias{logDiscr}
\alias{bifAttr}

\title{
Population Dynamic Models
}

\description{
Functions to simulate population dynamic models. 
}

\usage{
popExp(N0,lamb,tmax, intt)
estEnv(N0,lamb,varr,tmax, npop, ext) 
estDem(tmax,n ,m ,N0 , nsim, ciclo)
popLog(N0, r, K, tmax)
popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, rw, cl, tmax)
}

\arguments{
  \item{N0}{number of individuals at start time.
}
  \item{lamb}{finite rate of population growth.
}
  \item{intt}{interval time size.
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
  \item{rd}{discrete growth rate.
}
  \item{minrd}{minimum discrete growth rate.
}
  \item{maxrd}{maximum discrete growth rate.
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
popExp simulates an exponential population growth with environmental stochasticity. 

estDem simulates birth and death rates for stochastic (demographic) models.

Stochastic models uses lambda values taken of a normal distribution with mean lambda and variance varr.

popLog simulates a logistic growth for continuous and discrete models.

popStr simulates a structured population dynamics, with Leftkovitch matrices.

In popStr the number of patches in the simulated scene is defined by rw*cl.

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
