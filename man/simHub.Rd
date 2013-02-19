\Name{simHub}
\Alias{simHub1}
\Alias{simHub2}
\Alias{simHub3}

\Title{ 
Neutral Theory of Biogeography
}

\Description{
Simulates the Neutral Theory of Biogeography by Hubbell.
}

\Usage{
simHub1(S, j, D, cicles, anima=TRUE)
simHub2(S, j, D, cicles, m, anima=TRUE)
simHub3(Sm, jm, S, j, D, cicles, m, nu, anima=TRUE)
}

\Arguments{
  \item{S}{
    number of species in the community.
}
  \item{Sm}{
    number of species in the metacommunity.
}
  \item{j}{
    individuals per species in the community.
}
  \item{jm}{
    individuals per species in the metacommunity.
}
  \item{D}{
    number of dead per cicle.
}
  \item{cicles}{
    number of cicles in the simulation.
}
  \item{m}{
    colonization/immigration rate.
}
  \item{nu}{
    speciation rate.
}
  \item{anima}{
    logical; if TRUE, the simulation frames of the metacommunity are shown.
}
}

\Details{

}

\Value{

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
simHub1(S=10,j=10, D=1, cicles=5e3)
simHub2(j=2,cicles=2e4,m=0.1)
simHub3(Sm=200, jm=20, S= 10, j=100, D=1, cicles=1e4, m=0.01, nu=0.001, anima=TRUE)
}

\Keyword{simulation}
\Keyword{Neutral Theory}
