\name{meta.spac}
\alias{meta.spac}
\title{
Espatial Dependence Metapopulation Simulation
}
\description{
Simulate a metapopulation dynamics with espatial dependence
}
\usage{
meta.spac(tf, cl, ln, fi, pe, pc, canto = FALSE)
}
\arguments{
   \item{tf}{final time, numeric
}
  \item{cl}{number of coluns
}
  \item{ln}{number of rowns
}
  \item{fi}{initial fraction of patchs occupaded)
}
  \item{pc}{probability of colonization
}
  \item{pe}{probability of extinction
}
  \item{canto}{logical argument indicating if the initial occurence is/isn't cluster at the arena corner  
}
}
\value{
array with simulation data and a graphic device representing it
}
\references{
Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule rain, and the core-satellite hypothesis. The American Naturalist 138:768-776.
}
\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}
\examples{
metapop(tf=100,cl=10,ln=10,fi=0.5,pc=0.3,pe=0.15)
}
\keyword{metapopulation}
\keyword{simulation}
