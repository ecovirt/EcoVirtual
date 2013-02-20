\Name{arquip}

\Title{
Species Colonization and Species-Area Relationship in Arquipelagos
}

\Description{
Simulate species colonization in islands with different sizes, from the mainland with a certain species pool abundance and seed rain.
}

\Usage{
arquip(n.isl, ar.min, ar.max, S, seed.rain, abund, tmax=100, anima=TRUE)
}

\Arguments{
  \item{n.isl}{
    number of islands.
}
  \item{ar.min}{
     area of the smallest island.
}
  \item{ar.max}{
      area of the biggest island.
}
  \item{S}{
     number of species (species richness from mainland).
}
  \item{seed.rain}{
     seed rain. Number of seeds colonizing islands on each time.
}
  \item{abund}{
     abundance of each species in the seed rain.
}
 \item{tmax}{
     maximum time for the simulations.
}
 \item{anima}{
     logical; if TRUE, show simulation frames.
}
}

\Details{
The 'abund' argument can be a vector with the same length of the species richness, or a value representing equal abundance of each species (maximum eveness), or if the number varies between 0 and 1, the model of geometric species rank-abundance distribution is performed to calculate the abundances. The model is: abund*(1-abund)*((1:S)-1), where S is the number of species.
}

\Value{
arquip returns 3 graphs with the simulation results:
- The species-area relationship: number of species x island area at the end of the simulation. It also returns the coefficients c and z.
- Colonization rate curves: colonization (species/cicle) x species number for each island.
- Passive colonization: number of species x time for each island.
}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>}

\References{

}

\Examples{

arquip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=10, tmax=100, anima=TRUE)

arquip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=0.5, tmax=100, anima=TRUE)

}

\Keyword{simulation}
\Keyword{island biogeography}
