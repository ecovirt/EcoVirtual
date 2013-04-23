\name{metapopulation}
\alias{metapopulation}
\alias{metaPop}
\alias{metaCi}
\alias{metaEr}
\alias{metaCiEr}

\title{
Metapopulation Models
}

\description{
Simulate metapopulation dynamics with propagulus seed rain, internal colonization and rescue effect.
}

\usage{
metaPop(cl,rw, fi, pc, pe, tmax)
metaCi(cl, rw, fi, ci, pe, tmax)
metaEr(cl, rw, fi, pc, ce, tmax)
metaCiEr(cl, rw, fi, ci, ce, tmax)
}

\arguments{
  \item{cl}{
  number of columns to build the simulated landscape.
}
  \item{rw}{
  number of rowns to build the simulated landscape.
}
  \item{fi}{
  initial proportion of occupied patches.
}
  \item{pc}{
  probability of colonization.
}
  \item{pe}{
  probability of extinction.
}
  \item{ci}{
  colonization coefficient.
}
  \item{ce}{
   extiction coefficient.
}
  \item{tmax}{
    maximum simulation time.
}
}

\details{
'metaPop' is the seed rain metapopulation model, including only propagulus seed rain from a external pool (no extinction).

'metaCi' is the Internal Colonization model, where number of propagulus depends on number of occupied patches, there is no external pool.

'metaEr' is the Rescue Effect model, where extinction probability is negatively associated with number of occupied patches. 

'metaCiEr' includes both effects: Rescue Effect and Internal Colonization.

The number of patches in the simulated landscape is defined by rw*cl. 
}

\value{
Metapopulation functions return graphics with the simulation results.
These functions also return an invisible array with the simulation data.

} 

\references{
Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule rain, and the core-satellite hypothesis. The American Naturalist 138:768-776.

Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\examples{
metaPop(tmax=100,cl=10,rw=10,fi=0.5,pc=0.3,pe=0.15)
metaCi(tmax=100,cl=10,rw=10,fi=0.1,ci=1,pe=0.5)
metaEr(tmax=100, cl=10, rw=10, fi=0.2, pc=0.2, ce=0.15)
metaCiEr(tmax=100, cl=10, rw=10, fi=0.2, ci=0.2, ce=0.15)
}

\keyword{metapopulation}
\keyword{simulation}
