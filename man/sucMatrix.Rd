\Name{sucMatrix}

\Title{
Successional Stages Matrix
}

\Description{

}

\Usage{
sucMatrix(mat.trans, init.prop, rw, cl, tmax)
}

\Arguments{
  \item{mat.trans}{
    matrix of transitions probabilites.
}
  \item{init.prop}{
    initial stages proportions.
}
  \item{rw}{
    number of rows.
}
  \item{cl}{
    number of columns.
}
  \item{tmax}{
    maximum simulation time.
}
}

\Details{

}

\Value{

}

\References{

}

\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{

}

\Examples{
sucMatrix(mat.trans=matrix(data=c(0.5,0.5,0.5,0.5), nrow=2), init.prop=c(0.5,0.5),rw=20,cl=20, tmax=100)

}

\Keyword{simulation}
\Keyword{ }
