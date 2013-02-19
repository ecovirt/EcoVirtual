\name{metapopulation}
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
  number of columns.
}
  \item{rw}{
  number of rowns.
}
  \item{fi}{
  initial fraction of patchs occupaded.
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
     extiction coeficient.
}
  \item{tmax}{
    maximum simulation time.
}
}

\details{

}

\value{
return a graphic with the simulation results. 
} 

\references{
Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule rain, and the core-satellite hypothesis. The American Naturalist 138:768-776.
}

\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\examples{
metapop(tmax=100,cl=10,rw=10,fi=0.5,pc=0.3,pe=0.15)
metaCi(tmax=100,cl=10,rw=10,fi=0.1,ci=1,pe=0.5)
metaEr(tmax=100, cl=10, rw=10, fi=0.2, pc=0.2, ce=0.15)
metaCiEr(tf=100, cl=10, ln=10, fi=0.2, i=0.2, e=0.15)
}

\keyword{metapopulation}
\keyword{simulation}
