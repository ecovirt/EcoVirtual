\Name{rich}

\Title{
Number of Species
}

\Description{
Count the number of species (species richness) from a vector with a species list.
}

\Usage{
rich(x)
}

\Arguments{
  \item{x}{
    a vector with names.
}

}

\Details{
This function is used internally in the functions 'simHub1', simHub2', and 'simHub3'.
}

\Value{
returns the number of species (species richness).
}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\Examples{
lsp <- sample(LETTERS,50,replace=TRUE)
lsp
rich(lsp)
}

\Keyword{simulation}
\Keyword{Neutral Theory}
