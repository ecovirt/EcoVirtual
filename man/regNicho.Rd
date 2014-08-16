\name{regNicho}
\alias{regNicho}

\title{
Successional Niche Model				
}

\description{
Simulates the process of niche succession by successional stages in a community with 2 species (a superior and an inferior competitor), following the model of Pacala and Rees (1998).
}

\usage{
regNicho(tmax, rw, cl, c1, c2, ec, dst, er, sc, mx, rs)
}

\arguments{
  \item{tmax}{
    maximum simulation time.
}
  \item{rw}{
    number of rows for the simulated landscape.
}
  \item{cl}{
    number of columns for the simulated landscape.
}
  \item{c1}{
    colonization rate for the late successional species (superior competitor).
}
  \item{c2}{
    colonization rate for the early successional species (inferior competitor).
}
  \item{ec}{
    rate of competitive exclusion.
}
  \item{dst}{
    disturbance rate.
}
  \item{er}{
    inicial proportion of patches in early stage.
}
  \item{sc}{
    inicial proportion of patches in susceptible stage.
}
  \item{mx}{
    initial proportion of patches in mixed stage.
}
  \item{rs}{
    initial proportion of pathces in resistant stage.
}
}

\details{
There are five possible states of this model:
\itemize{
 \item \emph{free} - open, unoccupied space; 
 \item \emph{early} - occupied by only the early successional species; 
 \item \emph{susceptible} - occupied by only the late succesional species and susceptible to invasion by the early succesional species; 
 \item \emph{mixed} - occupied by both species; 
 \item \emph{resistant} - occupied by only the late successional species. 
}

The early successional species is the inferior competitor in the model, and the later successional species is the superior competitor.

The number of patches in the simulated landscape is defined by rw*cl. 

'dst' is the proportion of patches in any stage that turns empty, it represents a disturbance in the landscape.

'ec' is the probability of succeptible and mixed stages turns resistant stage.

}

\value{
'regNicho' returns the simulation results of patch occupancy in time for each successional stage.

'regNicho' also returns an invisible array with the simulation results per time.

}

\references{
Pacala, S & Rees, M. 1998. Models suggesting field experiments to test two hypotheses explaining successional diversity. The American Naturalist 152(2): 729:737.

Stevens, MHH. 2009. A primer in ecology with R. New York, Springer.
}

\author{
Alexandre Adalardo de Oliveira \email{ecovirtualpackage@gmail.com}
}

\seealso{
	\code{\link{comCompete}},
	 \url{http://ecovirtual.ib.usp.br}
}

\examples{
\dontrun{
regNicho(tmax=50, rw=100, cl=100, c1=0.2, c2=0.8, ec=0.5, dst=0.04,  er=0.08, sc=0.02, mx=0, rs=0)
}
}

\keyword{simulation}
\keyword{Niche sucession}
