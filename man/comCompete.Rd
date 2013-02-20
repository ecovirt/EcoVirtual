\Name{comCompete}

\Title{
Trade-off between colonization and competition species features
}

\Description{
Simulates the trade-off between colonization and competition features of species in the proportion of patch occupancy.
}

\Usage{
comCompete(rw, cl, S, fi, fsp1, pe, fr=0, int=0, tmax)
}

\Arguments{
  \item{rw}{
    number of rows.
}
  \item{cl}{
    number of columns.
} 
  \item{S}{
    number of species.
}
  \item{fi}{
    initial fraction of patchs occupaded.
}
  \item{fsp1}{
    superior competitor abundance.
}
  \item{pe}{
    mortality rate.
}
  \item{fr}{
    disturbance frequency.
}
  \item{int}{
    disturbance intensity.
}
  \item{tmax}{
    maximum simulation time.
}
}

\Details{

}

\Value{
return the graph of simulation results, with the frequency of patch occupied in time by each species, and the trade-off scale, with the superior competitor in one side and the superior colonizator in the other.
}

\References{

}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{

}

\Examples{
comCompete(tmax=1000, rw=100, cl=100, S=10, fi=1, fsp1=0.20, pe=0.01, fr=0, int=0)
}

\Keyword{simulation}
\Keyword{metacompetition}
