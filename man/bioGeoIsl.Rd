\name{bioGeoIsl}
\alias{bioGeoIsl}

\title{
Island Biogeographical Model			
}

\description{
Simulates island biogeographical models, with rates of colonization and extinction for islands of different sizes and distances to the mainland.
}

\usage{
bioGeoIsl(area, dist, P, weight.A=0.5, a.e=1, b.e=-0.01, c.i=1, d.i=-0.01, 
e.i=0, f.i=.01, g.e=0, h.e=.01)
}

\arguments{
  \item{area}{
    a vector with the sizes of the island areas. It must have the same length as 'dist'
}
  \item{dist}{
    a vector with the distances of the islands to the mainland. It must have the same length as 'areas'
}
  \item{P}{
    the number of species in the mainland (species richness of the pool).
}
  \item{weight.A}{
    ratio between the area and the distance effects.
}
  \item{a.e}{
    basal extiction coefficient for area.
}
  \item{b.e}{
   extinction/area coefficient.
}
  \item{c.i}{
    basal colonization coefficient for distance.
}
  \item{d.i}{
    numeric, colonization/distance coefficient.
}
  \item{e.i}{
    basal colonization coefficient for area.
}
  \item{f.i}{
   colonization/area coefficient.
}
  \item{g.e}{
    basal extinction coefficient for distance.
}
  \item{h.e}{
    extinction/distance coefficient.
}
}

\value{
'bioGeoIsl' returns a graph with the rates of colonization and extinction in relation with the species richness for each island.

'bioGeoIsl' also returns a invisible data frame with the values for area, distance and species richness (S) for each island.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\seealso{
	\code{\link{animaColExt}}
	\code{\link{archip}}
}


\examples{
\dontrun{
bioGeoIsl(area=c(5,10,50,80), dist=c(10,100,100,10), P=100, weight.A=.5, a=1, 
b=-0.01, c=1, d=-0.01, e=0, f=.01, g=0, h=.01)
}
}

\keyword{simulation}
\keyword{island biogeography}
