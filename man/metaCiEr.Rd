\name{metaCiEr}
\alias{metaCiEr}
\title{
Internal Colonization And Rescue Effect Simulation
}
\description{
Simulate a population with internal colonization and rescue effect
}
\usage{
metaCiEr(tf, cl, ln, fi, i, e)
}
\arguments{
  \item{tf}{
    simulation time
}
  \item{cl}{
     arena number of colunms
}
  \item{ln}{
      arena number of rowns
}
  \item{fi}{
     fraction of patches ocuppied
}
  \item{i}{
     colonization coeficient
}
  \item{e}{
     extinction coeficiente
}
}
\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>}
\examples{
metaCiEr(tf=100, cl=10, ln=10, fi=0.2, i=0.2, e=0.15)
}
\keyword{simulation}
\keyword{population}
