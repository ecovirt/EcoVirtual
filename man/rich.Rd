\name{rich}
\alias{rich}

\title{
Number of Species
}

\description{
Count the number of species (species richness) from a vector with a species list.
}

\usage{
rich(x)
}

\arguments{
  \item{x}{
    a vector with names.
}

}

\details{
This function is used internally in the functions 'simHub1', simHub2', and 'simHub3'.
}

\value{
returns the number of species (species richness).
}

\author{
Alexandre Adalardo de Oliveira \email{aleadalardo@gmail.com}
}

\examples{
lsp <- sample(LETTERS,50,replace=TRUE)
lsp
rich(lsp)
}

\keyword{simulation}
\keyword{Neutral Theory}
\keyword{Internal Functions}
