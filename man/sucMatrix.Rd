\name{sucMatrix}
\alias{sucMatrix}

\title{
Successional Stages Matrix
}

\description{
Simulates a successional model based on a transitional matrix of stages and its initial proportion of occurence in the landscape.
}

\usage{
sucMatrix(mat.trans, init.prop, rw, cl, tmax)
}

\arguments{
  \item{mat.trans}{
   a matrix of stage transition probabilites.
}
  \item{init.prop}{
    a vector with the initial proportions of each stage.
}
  \item{rw}{
    number of rows to build the simulated landscape.
}
  \item{cl}{
    number of columns to build the simulated landscape.
}
  \item{tmax}{
    maximum simulation time.
}
}

\details{


The number of patches in the simulated landscape is defined by rw*cl.  
}

\value{
'sucMatrix' return a simulation graphic with the proportions of stages in the landscape in time, and a stage distribution graphic with the results of the simulation with the number o patches in time for each stage.

'sucMatrix' also return an invisible array with the simulation results.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}


\examples{
\dontrun{
sucMatrix(mat.trans=matrix(data=c(0.5,0.5,0.5,0.5), nrow=2), 
init.prop=c(0.5,0.5),rw=20,cl=20, tmax=100)
}

}

\keyword{simulation}
\keyword{ecological succession} 
