\name{anima}
\alias{animaCena}
\alias{animaGame}
\alias{animaHub}
\alias{animaIsl}
\alias{animaMeta2}
\alias{animaMetaComp}
\alias{animaRandWalk}
\alias{grColExt}
\alias{grFim}

\title{
EcoVirtual Graphics and Animations
}

\description{
Graphics and animations of the simulations results.
}

\usage{
animaCena(dados)
animagGame(xGame, total, sleep=0.01)
animaHub(dadoHub, sleep=0.1)
animaIsl(riq.tempo, ar.isl, locxy, sprain, col_riq=col_riq)
animaMeta2(dados)
animaMetaComp(dados)
animaRandWalk(rwData, time=2, sleep=0.1)
grColExt(E , I , P, areas)
grFim(dados)

}

\arguments{
  \item{dados}{

}
  \item{xGame}{

}
  \item{total}{

}
  \item{sleep}{

}
 \item{dadoHub}{

}
 \item{riq.tempo}{

}
 \item{ar.isl}{

}
 \item{locxy}{

}
 \item{sprain}{

}
 \item{col_riq}{

}
 \item{rwdata}{

}
 \item{time}{

}
}

\details{
Below each function and its primary function (where the graphical and animation function are used):
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
grColExt(E = 0.5 , I = 0.5 , P = 100, areas=1:10)
}
\keyword{ simulation }
