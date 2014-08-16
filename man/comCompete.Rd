\name{comCompete}
\alias{comCompete}

\title{
Multispecies competition-colonization tradeoff
}

\description{
Simulates the trade-off between colonization and competition abilities in a multispecies system.
}

\usage{
comCompete(rw, cl, S, fi, fsp1, pe, fr=0, int=0, tmax)
}

\arguments{
  \item{rw}{
    number of rows for the simulated landscape.
}
  \item{cl}{
   number of columns for the simulated landscape.
} 
  \item{S}{
    number of species.
}
  \item{fi}{
    initial fraction of patchs occupied
}
  \item{fsp1}{
    superior competitor abundance.
}
  \item{pe}{
    mortality rate.
}
  \item{fr}{
    disturbance frequency.
}
  \item{int}{
    disturbance intensity.
}
  \item{tmax}{
    maximum simulation time.
}
}

\details{
In the system, the competitive abilities are inversely proportional to the colonization abilities.

The number of patches in the simulated landscape is defined by rw*cl. 
}

\value{
'comCompete' returns a graph with the proportion of patches occupied in time by each species and the trade-off scale, the superior competitor in one side and the superior colonizator in the other.
}

\references{
Tilman. R. 1994. Competition and biodiversity in spatially structured habitats. Ecology,75:2-16. 

Stevens, M.H.H. 2009. A primer in ecology with R. New York, Springer.
}

\author{
Alexandre Adalardo de Oliveira \email{ecovirtualpackage@gmail.com}
}

\seealso{
	\code{\link{metaComp}},
        \url{http://ecovirtual.ib.usp.br}

}

\examples{
\dontrun{
comCompete(tmax=1000, rw=100, cl=100, S=10, fi=1, fsp1=0.20, pe=0.01, fr=0, int=0)
}
}

\keyword{simulation}
\keyword{metacompetition}
