\name{randWalk}
\alias{randWalk}

\title{
Random Walk Simulations
}

\description{
Simulates random walk models.
}

\usage{
randWalk(S,step=1,tmax=1e5,x1max=200, alleq=FALSE)
}

\arguments{
  \item{S}{
    number of individuals.
}
  \item{step}{
    step size (number of steps on each time)
}
  \item{tmax}{
    maximum simulation time.
}
  \item{x1max}{
    maximum initial distance from absorption surface.
}
  \item{alleq}{
    logical; if TRUE, all initial distance are equal. if FALSE, initial distances for each individual is a sample between 1 and maximum initial distance(x1max).
}
}

\details{
Random walk is a stochastic process of a succession of random steps. 

Zero is the absorption surface. When an individual simulation reaches zero, it means that the individual is dead.

See \url{ http://en.wikipedia.org/wiki/Random_walk}.
}

\value{
'randWalk' returns a graphic with the simulated trajectories of each individual. 

'randWalk' also returns an invisible matrix with the distance from de edge for each individual on each time.
}

\references{
\url{ http://en.wikipedia.org/wiki/Random_walk}
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\seealso{
	\code{\link{extGame}}
  	\code{\link{simHub}}
}

\examples{
\dontrun{
randWalk(S=100,step=2,tmax=2e5)
randWalk(S=10,step=1,tmax=1e4, x1max=300, alleq=TRUE)
}
}

\keyword{simulation}
\keyword{Neutral Theory}
