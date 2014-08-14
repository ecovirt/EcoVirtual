\name{archip}
\alias{archip}

\title{
Species Colonization and Species-Area Relationship in Archipelagos
}

\description{ 
Simulate species colonization from mainland to islands with different sizes. 
}

\usage{
archip(n.isl, ar.min, ar.max, S, seed.rain, abund, tmax=100, anima=TRUE)
}

\arguments{
  \item{n.isl}{
    numeric, number of islands.
}
  \item{ar.min}{
     numeric, area of the smallest island.
}
  \item{ar.max}{
     numeric, area of the biggest island.
}
  \item{S}{
     numeric, number of species (species richness from mainland).
}
  \item{seed.rain}{
    numeric, seed rain. Number of seeds colonizing islands on each time.
}
  \item{abund}{
    numeric, abundance of each species in the seed rain.
}
 \item{tmax}{
    numeric, maximum time for the simulations.
}
 \item{anima}{
     logical; if TRUE, show simulation frames.
}
}

\details{
The mainland has richness (S) and the evenness can be controled argument abund. The 'abund' argument can be one of these 3 options:
\enumerate{
	\item a vector with the same length of the species richness, meaning the proportion of each species population; 
	\item a single value more than 1 representing equal abundance of each species (maximum evenness); 
	\item a single value between 0 and 1, meaning the model of geometric species rank-abundance distribution. The model is: abund*(1-abund)*((1:S)-1), where S is the number of species.
}
}

\value{
'archip' returns 3 graphics:
\itemize{
	\item The species-area relationship: number of species x island area at the end of the simulation. It also returns the coefficients c and z from species-area relationship \eqn{S=cA^z}.

	\item Colonization rate curves: colonization (number of species per cycle) x number of species for each island.

	\item Passive colonization: number of species x time for each island.

'archip' also returns an invisible array with the simulation results.
}
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\seealso{
	\code{\link{animaColExt}},
	\code{\link{bioGeoIsl}}
}

\examples{

\dontrun{
archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=10, tmax=100, anima=TRUE)

archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=0.5, tmax=100, anima=TRUE)
}

}

\keyword{simulation}
\keyword{island biogeography}
\keyword{Species-area relationship}
