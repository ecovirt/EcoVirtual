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
estExp(N0,r,varr,tmax) 
estDem(N0, b, d, tmax)
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
  \item{p.aa}{adult probability of survival.
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
'estLog' give results for continuous and discrete models.
'estExp' and 'estDem' give results for deterministic and stochastic models.

}

}
\value{
return a grafic with the simulation results.
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
