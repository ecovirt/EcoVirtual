\name{animaColExt}
\alias{animaColExt}

\title{ 
Colonization and Extinction balance in the Island Biogeography Equilibrium model
}

\description{
Simulate the balance between extinction and colonization rates given the equilibrium number of species in a island, based on the Island Biogeography Equilibrium model.
}

\usage{
animaColExt(min=0.01, max=1, cycles=100, Ext="crs", Col="dcr")
}

\arguments{
  \item{min}{ 
    between 0-1. The minimum value of the extinction and colonization rates.
}
  \item{max}{
    between 0-1. The maximum value of the extinction and colonization rates.
}
  \item{cycles}{
    number of cycles in the simulation.
}
  \item{Ext}{
   a string representing the extinction rate. This can be 'crs' for an increasing extinction rate, 'fix' for a fixed extinction rate in 0.5, or 'dcr' for a decreasing extinction rate.   
}
  \item{Col}{
   a string representing the colonization rate. This can be 'crs' for an increasing colonization rate, 'fix' for a fixed colonization rate in 0.5, or 'dcr' for a decreasing colonization rate.
}
}

\details{
 The number of species is the balance between extinction and colonization rates at the equilibrium.		
}

\value{
'animaColExt' returns a graph of the extinction and colonization rates varying one or both rates in relation with the number of species of an island.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}}

\seealso{
	\code{\link{arquip}}
	\code{\link{bioGeoIsl}} 
}

\examples{
\dontrun{
animaColExt(Ext='fix', Col="fix")
}
}

\keyword{simulation}
\keyword{Island Biogeography}

