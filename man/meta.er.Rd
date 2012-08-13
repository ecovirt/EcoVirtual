\name{meta.er}
\alias{meta.er}
\title{
  Rescue Effect Simulation
}
\description{
  Simulate a population with rescue effect
}
\usage{
meta.er(tf, cl, ln, fi, pc, e)
}
\arguments{
  \item{tf}{
    simulation time
}
  \item{cl}{
     arena number of colunms
}
  \item{ln}{
      arena number of rowns
}
  \item{fi}{
     fraction of patches ocuppied
}
  \item{pc}{
     colonization probability
}
  \item{e}{
		extiction coeficient
}
}
\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>}
\examples{
meta.cier(tf=100, cl=10, ln=10, fi=0.2, i=0.2, e=0.15)
}
\keyword{simulation}
\keyword{population}
