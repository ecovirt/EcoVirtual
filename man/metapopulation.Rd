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
Simulate metapopulation dynamics with propagules seed rain, internal colonization and rescue effect.
}

\usage{
metaPop(cl,rw, f0, pi, pe, tmax)
metaCi(cl, rw, f0, ci, pe, tmax)
metaEr(cl, rw, f0, pi, ce, tmax)
metaCiEr(cl, rw, f0, ci, ce, tmax)
}

\arguments{
  \item{cl}{
  number of columns for the simulated landscape.
}
  \item{rw}{
  number of rows for the simulated landscape.
}
  \item{f0}{
  initial proportion of occupied patches.
}
  \item{pi}{
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
'metaPop' is the seed rain metapopulation model, including only propagules seed rain from a external pool (no extinction).

'metaCi' is the Internal Colonization model, where number of propagules depends on number of occupied patches, there is no external pool.

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
Alexandre Adalardo de Oliveira and Paulo Inacio Prado \email{ecovirtualpackage@gmail.com}
}

\seealso{
        \url{http://ecovirtual.ib.usp.br}
}
\examples{
\dontrun{
metaPop(cl=10,rw=10,f0=0.5,pi=0.3,pe=0.15, tmax=100)
metaCi(cl=10,rw=10,f0=0.1,ci=1,pe=0.5, tmax=100)
metaEr(cl=10, rw=10, f0=0.2, pi=0.2, ce=0.15, tmax=100)
metaCiEr(cl=10, rw=10, f0=0.2, ci=0.2, ce=0.15, tmax=100)
}
}

\keyword{metapopulation}
\keyword{simulation}
