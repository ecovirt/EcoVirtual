\name{metaInter}
\alias{metaInter}
\title{
Internal Colonization 
}
\description{
Model a population with internal colonization
}
\usage{
metaInter(tf, cl, ln, fi, i, pe)
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
  \item{i}{colonization constante 
}
  \item{pe}{probability of extinction
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
metaInter(tf=100,cl=10,ln=10,fi=.1,i=1,pe=0.5)
}
\keyword{metapopulation}
\keyword{simulation}
