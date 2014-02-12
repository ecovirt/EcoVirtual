\name{dynPop}
\alias{popExp}
\alias{estEnv}
\alias{BDM}
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
popExp(N0, lamb, tmax, intt)
estEnv(N0, lamb, tmax, varr, npop, ext) 
BDM(tmax, b, d, m, n0)
estDem(N0, tmax, n ,m , migr, nsim, ciclo)
popLog(N0, tmax, r, K, ext)
popStr(tmax, p.sj, p.jj, p.ja, p.aa, fec, ns, nj, na, rw, cl)
logDiscr(N0, rd, K, tmax, ...)
}

\arguments{
  \item{N0}{number of individuals at start time.
}
  \item{lamb}{finite rate of population growth.
}
  \item{tmax}{maximum simulation time.
}
  \item{intt}{interval time size.
}
  \item{varr}{variance.
}
  \item{npop}{
}
  \item{ext}{
}
  \item{b}{
}
  \item{d}{
}
  \item{n0}{
}
  \item{n}{birth rate.
}
  \item{m}{death rate.
}
  \item{migr}{
}
  \item{nsim}{number of populations simulated.
}
  \item{ciclo}{size of cicles in simulation.
}
  \item{r}{intrinsic growth rate.
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
  \item{fec}{mean number of propagulus per adult each cicle.
}
  \item{ns}{number of seeds at initial time.
}
  \item{nj}{number of juvenil at initial time.
}
  \item{na}{number of adults at initial time.
}
  \item{rw}{number of rows to build the simulated scene.
}
  \item{cl}{number of columns to build the simulated scene.
}
  \item{rd}{
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
