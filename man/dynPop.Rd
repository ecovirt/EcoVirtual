\name{dynPop}
\alias{estExp}
\alias{estDem}
\alias{crescLog}
\alias{popStr}

\title{
Populations Dynamics Model
}

\description{
Functions to simulate population dynamics 
}

\usage{
estExp(N0,r,varr,tmax) 
estDem(N0, b, d, tmax)
crescLog(N0, r, K, tmax)
popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, ln, cl, tmax)
}

\arguments{
  \item{N0}{number of individuos at start time
}
  \item{tmax}{final time, numeric
}
  \item{cl}{number of coluns
}
  \item{ln}{number of rowns
}
  \item{r}{intrinsic growth rate)
}
  \item{varr}{r variance
}
  \item{b}{birth rate
}
  \item{d}{dead rate
}
  \item{K}{carrying capacity
}
  \item{p.sj}{probability of seed survivel
}
  \item{p.jj}{probability of juvenil survivel
}
  \item{p.ja}{juvenil probability of transition to adult fase
}
  \item{p.aa}{adult probability of survivel 
}
  \item{ns}{number of seeds at initial time
}
  \item{na}{number of adults at initial time
}
  \item{nj}{number of juvenil at initial time
}
  \item{fec}{mean number of propagulus per adult each cicle
}

}
\value{
return a grafic with the simulation and an array with all simulation data
}
\references{
Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule rain, and the core-satellite hypothesis. The American Naturalist 138:768-776.
}
\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}
\examples{
popStr(p.sj=0.4, p.jj=0.6, p.ja=0.2, p.aa=0.9, fec=0.8, ns=100,nj=40,na=20, ln=30, cl=30, tmax=20)
}
\keyword{population dynamics}
\keyword{simulation}
