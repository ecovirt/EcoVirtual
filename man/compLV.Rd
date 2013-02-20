\Name{compLV}

\Title{
Lotka-Volterra Competition Model	
}

\Description{
Simulate the Lotka-Volterra competition model for two populations.
}

\Usage{
compLV(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
}

\Arguments{
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

\Details{

}

\Value{
'compLV' returns a graph of the population size in time, and a graph with the isoclines of the equilibrium for both species.
}

\References{
Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer Associates, 291pp.
}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{

}

\Examples{
compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)
}

\Keyword{simulation}
\Keyword{competition}
