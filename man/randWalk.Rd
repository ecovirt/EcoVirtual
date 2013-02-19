\Name{randWalk}

\Title{
Random Walk Simulation
}

\Description{
Simulate the extinction of populations by a randow walk model.
}

\Usage{
randWalk(S,step=1,tmax=1e5,x1max=200, all.eq=FALSE)
}

\Arguments{
  \item{S}{
    number of species or populations.
}
  \item{step}{
    step size. number of death or birth on each population in each time.
}
  \item{tmax}{
    maximum time.
}
  \item{x1max}{
    maximum initial distance to extinction for each population.
}
  \item{all.eq}{
    logical; if TRUE, the initial distance of the species are equal. if FALSE, the initial distance of the species are a sample between 1 and the maximum initial distance(x1max).
}
}

\Details {
Step size determines the amount of birth or death on each time, by chance. The distance from the edge tels the population.
}

\Value{
randWalk returns a graph with the trajectories of the populations in time. Extincted populations are those who reached the absortion surface.
}

\References{

}

}
\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{

}

\Examples{
randwalk(S=100,step=2,tmax=2e5)
randWalk(S=10,step=1,tmax=1e4, x1max=300, all.eq=TRUE)

}

\Keyword{simulation}
\Keyword{Neutral Theory}
