\name{extGame}
\alias{extGame}

\title{
Zero-sum game
}

\description{
Simulates a zero-sum game between two competitors with a fixed amount of resource.
}

\usage{
extGame(bet, total, tmax)
}

\arguments{
  \item{bet}{
    bet size of each competitor on each time.
}
  \item{total}{
    total amount of resource.
}
  \item{tmax}{
    maximum game time.
}
}

\details{
A zero-sum game is a mathematical representation of a situation in which a participant's gain (or loss) of resource is exactly balanced by the losses (or gains) of the resource of the other participant(s). If the total gains of the participants are added up, and the total losses are subtracted, they will sum to zero. 
}

\value{
'extGame' returns a graphic with the amount of resource of each competitor on each time.

'extGame' also returns an invisible vector with the results of the loser on each time.
}

\references{
\url{http://en.wikipedia.org/wiki/Zero-sum_game}
}

\author{
Alexandre Adalardo de Oliveira and Paulo Inacio Prado \email{ecovirtualpackage@gmail.com}
}


\seealso{
	\code{\link{simHub}},
  	\code{\link{randWalk}},
	\url{http://ecovirtual.ib.usp.br}
}

\examples{
\dontrun{
extGame(bet=1,total=20)
extGame(bet=1,total=100)
}
}

\keyword{simulation}
\keyword{neutral theory}
