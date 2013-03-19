\Name{randWalk}

\Title{
Random Walk Simulations
}

\Description{
Simulates random walk models.
}

\Usage{
randWalk(S,step=1,tmax=1e5,x1max=200, all.eq=FALSE)
}

\Arguments{
  \item{S}{
    number of individuals.
}
  \item{step}{
    step size. number of steps on each time.
}
  \item{tmax}{
    maximum time.
}
  \item{x1max}{
    maximum initial distance to each individual.
}
  \item{all.eq}{
    logical; if TRUE, the all initial distance are equal. if FALSE, the initial distances for each individual is a sample between 1 and the maximum initial distance(x1max).
}
}

\Details {
Random walk is a mathematical formalization of a path that consists of a succession of random steps. 
}

\Value{
'randWalk' returns a graph with the trajectories of the populations in time. Extincted populations are those who reached the absortion surface, when population size is zero.

'randWalk' also returns an invisible matrix with the distance from de edge for each individual on each time.
}

\References{

}

}
\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{
extGame, simHub
}

\Examples{
randwalk(S=100,step=2,tmax=2e5)
randWalk(S=10,step=1,tmax=1e4, x1max=300, all.eq=TRUE)

}

\Keyword{simulation}
\Keyword{Neutral Theory}
