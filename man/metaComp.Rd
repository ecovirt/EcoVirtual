\name{metaComp}

\title{
Metapopulation Competition Model
}

\description{
Simulate a metapopulation dynamics with two species with differente competition ability
}

\usage{
metaComp(tmax, rw, cl, fi1, fi2, i1, i2, pe, D=0, anima=TRUE)
}

\arguments{
  \item{tmax}{
maximum simulation time.
}
  \item{cl}{
number of columns.
}
  \item{rw}{
number of rowns.
}
  \item{fi1}{
initial fraction of patchs occupeded by sp1.
}
  \item{fi2}{
initial fraction of patchs occupeded by sp2.
}
  \item{i1}{
colonization coefficient for sp1.
}
  \item{i2}{
colonization coefficient for sp2.
}
  \item{pe}{
probability of extinction.
}
  \item{D}{
habitat destruction, proportion.
}
  \item{anima}{
show animation,  logical.
}
}

\details{

}

\value{
return a grafic with the simulation
}

\references{
Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule rain, and the core-satellite hypothesis. The American Naturalist 138:768-776.
}
\author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\examples{
metaComp(tmax=100,cl=20,rw=20,fi1=0.1,fi2=0.4,i1=0.4,i2=0.5,pe=0.25)
}

\keyword{metapopulation}
\keyword{simulation}
