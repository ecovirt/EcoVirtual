\name{anima}
\alias{animaCena}
\alias{animaGame}
\alias{animaHub}
\alias{animaIsl}
\alias{animaMeta2}
\alias{animaMetaComp}
\alias{animaRandWalk}
\alias{grColExt}
\alias{grToff}
\alias{grFim}

\title{
EcoVirtual Graphics and Amimations
}

\description{
Graphics and animations of the simulations results
}

\usage{
grToff(rq, fsp1,pe,add=FALSE, ...)
}

\arguments{
  \item{rq}{
number of species
}
  \item{fsp1}{
initial proportion for the most abundance species 
}
  \item{pe}{
extinction probability
}
  \item{add}{plot in the actual device or in a new one, logical 
}
 \item{...}{
another argments passed to grToff
}
}

\value{
Represent the simulation result in a graphic device
}

\author{
Alexandre Adalardo
}

\examples{
grToff(rq =  10 , fsp1 =  0.2 , pe =  0.1 ,add=FALSE)
}
\keyword{ simulation }
