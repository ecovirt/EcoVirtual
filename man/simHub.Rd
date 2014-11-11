\name{simHub}
\alias{simHub}
\alias{simHub1}
\alias{simHub2}
\alias{simHub3}

\title{ 
Neutral Theory of Biogeography
}

\description{
Simulates Community Dynamics as in the Neutral Theory of Biogeography
}

\usage{
simHub1(S, j, D, cycles, m.weights=1, anima=TRUE)
simHub2(S, j, D, cycles, m, anima=TRUE)
simHub3(Sm, jm, S, j, D, cycles, m, nu, anima=TRUE)
}

\arguments{
  \item{S}{
    number of species in the community.
}
  \item{Sm}{
    number of species in the metacommunity.
}
  \item{j}{
    individuals per species in the metacommunity.
}
  \item{jm}{
    individuals per species in the metacommunity.
}
  \item{D}{
    number of deaths per cycle.
}
  \item{cycles}{
    number of cycles in the simulation.
}
  \item{m}{
    colonization/immigration rate.
}
  \item{nu}{
    speciation rate.
  }
  \item{m.weights}{
    Mortality weights for each species. Mortality rates of individuals of each species is
    proportional to species' abundances multiplied by these weights as
    in Yu et al. (1998). In neutral dynamics all weigths are equal. If
    length(m.weights)<S then species are divided in groups of
    (approximately) S/length(m.weights) and species of each group have a
    value in m.weights. This allows to create groups of species with
  different mortality probabilities and compare to the neutral dynamics.
    }
  \item{anima}{
    logical; if TRUE, the simulation frames of the metacommunity are shown.
}
}

\details{
'simHub1' is the model without immigration.

'simHub2' incorporates immigration rate from the metacommunity

'simHub3' incorporates immigration and speciation rates in the metacommunity.
}

\value{
These functions returns a graph with the number of species in time (cycles) in the metacommunity.

They also return an invisible matrix with the results of species richness on each community per time.
}

\references{
Hubbell, S.P. 2001. The Unified Neutral Theory of Biodiversity and
Biogeography. Princeton University Pres, 448p.

Yu, D. W., Terborgh, J. W., and Potts, M. D. 1998. Can high tree species
richness be explained by Hubbell's null model?. Ecology Letters, 1(3):
193--199. 
}


\author{
Alexandre Adalardo de Oliveira and Paulo Inacio Prado \email{ecovirtualpackage@gmail.com}
}

\seealso{
	\code{\link{extGame}},
  	\code{\link{randWalk}},
        \url{http://ecovirtual.ib.usp.br}
}

\examples{
\dontrun{
simHub1(S=10,j=10, D=1, cycles=5e3)
simHub2(j=2,cycles=2e4,m=0.1)
simHub3(Sm=200, jm=20, S= 10, j=100, D=1, cycles=1e4, m=0.01, nu=0.001, anima=TRUE)
}
}

\keyword{simulation}
\keyword{Neutral Theory}
