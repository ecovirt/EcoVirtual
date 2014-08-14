\name{anima}
\alias{animaCena}
\alias{animaGame}
\alias{animaHub}
\alias{animaIsl}
\alias{animaMeta2}
\alias{animaMetaComp}
\alias{animaRandWalk}
\alias{grColExt}
\alias{gr.toff}
\alias{grFim}
\alias{anima}

\title{
Internal EcoVirtual Graphics and Animations
}

\description{
Internal functions for graphics and animations of the simulations results. 
}

\usage{
grColExt(E , I , P, area)
}

\arguments{
  \item{E}{
	extinction rate
}
  \item{I}{
	colonization rate
}
  \item{P}{

}
  \item{area}{

}

}


\details{

	The list below relates each function graphical and its primary functions:

animaCena - regNicho

animaGame - extGame

animaHub - simHub1, simHub2, simHub3

animaIsl - arquip

animaMeta2 - metaPop, metaCi, metaEr, metaCiEr

animaMetaComp - metaComp

animaRandWalk - randWalk

grColExt - animaColExt, bioGeoIsl

grFim - metaPop, metaCi, metaEr, metaCiEr
}

\value{
Represent the simulation result in a graphic device.
}

\author{
Alexandre Adalardo
}

\examples{
\dontrun{
grColExt(E = 0.5 , I = 0.5 , P = 100, area=1:10)
}
}
\keyword{ simulation }
