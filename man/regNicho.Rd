\Name{regNicho}

\Title{
Niche Regeneration Model				
}

\Description{
Simulates the process of niche succession in a landscape, with the proportion of each succecional stage in time.
}

\Usage{
regNicho(tmax, rw, cl, c1, c2, ec, dst, er, sc, mx, rs)
}

\Arguments{
  \item{tmax}{
    maximum simulation time.
}
  \item{rw}{
    number of rows.
}
  \item{cl}{
    number of columns.
}
  \item{c1}{
    colonization rate for the superior competitor.
}
  \item{c2}{
    colonization rate for the inferior competitor.
}
  \item{ec}{
    competitive exclusion.
}
  \item{dst}{
    disturbance.
}
  \item{er}{
    early stage (only sp2).
}
  \item{sc}{
    susceptive stage (only sp1).
}
  \item{mx}{
    mixed (sp1 and sp2).
}
  \item{rs}{
    resistant stage.
}
}

\Details{

}

\Value{
returns the simulation results of patch occupancy in time for each succecional stage.
}

\References{
Pacala, S & Rees, M. 1998. Models suggesting field experiments to test two hypotheses explaining successional diversity. The American Naturalist 152(2): 729:737.
Stevens, MHH. 2009. A primer in ecology with R. New York, Springer.
}

}
\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\Seealso{
	\code{\link{comCompete}}
}

\Examples{
regNicho(tmax=50, rw=100, cl=100, c1=0.2, c2=0.8, ec=0.5, dst=0.04,  er=0.08, sc=0.02, mx=0, rs=0)
}

\Keyword{simulation}
\Keyword{Niche sucession}
