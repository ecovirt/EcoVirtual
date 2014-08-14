\name{compLV}
\alias{compLV}

\title{
Lotka-Volterra Competition Model	
}

\description{
Simulate the Lotka-Volterra competition model for two populations.
}

\usage{
compLV(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
}

\arguments{
  \item{n01}{
    initial population for the superior competitor species.
}
  \item{n02}{
    initial population for the inferior competitor species.
}
  \item{tmax}{
    maximum simulation time.
}
  \item{r1}{
    intrinsic growth rate for the superior competitor species.
}
  \item{r2}{
    intrinsic growth rate for the inferior competitor species.
}
  \item{k1}{
    carrying capacity for the superior competitor species.
}
  \item{k2}{
    carrying capacity for the inferior competitor species.
}
  \item{alfa}{
    alfa coefficient.
}
  \item{beta}{
    beta coefficient
}
}

%\details{
%The Lotka-Volterra competition model follows the equations:
%\deqn{	\item[-] SP1:
%	$$\frac{dN_1}{dt}=r_1N_1\left(\frac{K_1-N_1-\alpha N_2}{K_1}\right) $$
%	\item[-] SP2:
%	$$\frac{dN_2}{dt}=r_2N_2\left(\frac{K_2-N_2-\beta N_1}{K_2}\right) $$
%}}
\details{
The Lotka-Volterra competition model follows the equations:
\itemize{	\item SP1:
dN_1/dt = r_1*N_1*((K_1-N_1-\alpha N_2)/K_1)
	\item SP2:
dN_2/dt = r_2*N_2*((K_2-N_2-\beta N_1)/K_2)

\value{
'compLV' returns a graph of the population size in time, and a graph with the isoclines of the equilibrium for both species.
'compLV' also returns an invisible matrix with the population size of each species in time.
}

\references{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
Hastings, A. 1980. Disturbance, coexistence, history and competition for space. Theoretical Population Biology, 18:363-373.
Stevens, M.H.H. 2009. A primer in ecology with R. New York, Springer.

}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\examples{
\dontrun{
	compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)
}
}

\keyword{simulation}
\keyword{competition}
