\Name{arquip}
\Alias{arquip}

\Title{
Species Colonization and Species-Area Relationship in Arquipelagoes
}

\Description{
Simulate species colonization in islands of different sizes, from the mainland with a certain species pool abundance and seed rain.
}

\Usage{
arquip(n.isl, ar.min, ar.max, S, seed.rain, abund, tmax, anima=TRUE)
}

\Arguments{
  \item{n.isl}{
    number of islands
}
  \item{ar.min}{
     area of the smallest island
}
  \item{ar.max}{
      area of the biggest island
}
  \item{S}{
     number of species (species richness from mainland)
}
  \item{seed.rain}{
     seed rain. Number of seeds colonizing islands on each time(cicle).
}
  \item{abund}{
     abundance of each species. AQUI TEM QUE EXPLICAR DIREITO.
}
 \item{tmax}{
     maximum time for the simulations. Is the same as number of cicles.
}
 \item{anima}{
     logical; if TRUE, show simulation frames; if FALSE hide simulation frames.
}
}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>}

\Examples{
arquip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=rep(10,1000), tmax=100, anima=TRUE)
}

\Keyword{simulation}
\Keyword{island biogeography}
