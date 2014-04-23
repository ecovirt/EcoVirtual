\name{simHub}
\alias{simHub}
\alias{simHub1}
\alias{simHub2}
\alias{simHub3}

\title{ 
Neutral Theory of Biogeography
}

\description{
Simulates the Neutral Theory of Biogeography by Hubbell.
}

\usage{
simHub1(S, j, D, cicles, anima=TRUE)
simHub2(S, j, D, cicles, m, anima=TRUE)
simHub3(Sm, jm, S, j, D, cicles, m, nu, anima=TRUE)
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
    number of dead per cicle.
}
  \item{cicles}{
    number of cicles in the simulation.
}
  \item{m}{
    colonization/immigration rate.
}
  \item{nu}{
    speciation rate.
}
  \item{anima}{
    logical; if TRUE, the simulation frames of the metacommunity are shown.
}
}

\details{
'simHub1' is the model without immigration in the metacommunity

'simHub2' incorporates the immigration rate in the metacommunity

'simHub3' incorporates the immigration and speciation rates in the metacommunity.
}

\value{
These functions returns a graph with the number of species in time (cicles) in the metacommunity.

They also return an invisible matrix with the results of species richness on each community per time.
}

\references{
Hubbell, S.P. 2001. The Unified Neutral Theory of Biodiversity and Biogeography. Princeton University Pres, 448p.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\seealso{
 	\code{\link{extGame}}
  	\code{\link{randWalk}}
}

\examples{
simHub1(S=10,j=10, D=1, cicles=5e3)
simHub2(j=2,cicles=2e4,m=0.1)
simHub3(Sm=200, jm=20, S= 10, j=100, D=1, cicles=1e4, m=0.01, nu=0.001, anima=TRUE)
}

\keyword{simulation}
\keyword{Neutral Theory}
