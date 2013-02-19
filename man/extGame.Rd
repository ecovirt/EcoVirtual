\Name{extGame}

\Title{
Zero-sum Game			
}

\Description{
Simulate a zero-sum game between two 'competitor species' in a community with a finite number of individuals.
}

\Usage{
extGame(bet,total, tmax)
}

\Arguments{
  \item{bet}{
    bet size. number of individuals 'converted' on the winner species per time (cicle).
}
  \item{total}{
    total amount. maximum number of individuals in the community.
}
} \item{tmax}{
    maximum game time (number of cicles).
}
}

\Details{

}

\Value{
extGame returns a graph with the number of individuals of each competitor on each cicle
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
extGame(bet=1,total=20)
extGame(bet=1,total=100)

}

\Keyword{simulation}
\Keyword{neutral theory}
