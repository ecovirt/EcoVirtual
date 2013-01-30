\name{metapop}
\alias{metapop}
\title{
Metapopulation Model
}
\description{
Simulate a metapopulation dynamics
}
\usage{
metapop(tf, cl, ln, fi, pc, pe)
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
metapop(tf=100,cl=10,ln=10,fi=0.5,pc=0.3,pe=0.15)
}
\keyword{metapopulation}
\keyword{simulation}
