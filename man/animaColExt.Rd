\Name{animaColExt}
\Alias{animaColExt}

\Title{ 
Colonization and Extinction balance in Island Biogeography
}

\Description{
Simulate the relationship between extinction and colonization rates for the species richness in Island Biogeography models.
}

\Usage{
animaColExt(min=0.01, max=1, cicles=100, Ext="crs", Col="dcr")
}

\Arguments{
  \item{min}{ 
    the minimum value of the extinction and colonization rates.
}
  \item{max}{
    the maximum value of the extinction and colonization rates.
}
  \item{cicles}{
    number of cicles in the simulation
}
  \item{Ext}{
   extinction rate. It can be 'crs' for an increasing extinction rate, 'fix' for a fixed extinction rate in 0.5, or 'dcr' for a decreasing extinction rate.
    
}
  \item{Col}{
   colonization rate. It can be 'crs' for an increasing colonization rate, 'fix' for a fixed colonization rate in 0.5, or 'dcr' for a decreasing colonization rate.
}
}

\Details{

}

\Value{
animaColExt returns a graph of the extinction and colonization rates on each cicle in relation with te number of species of an island in Island Biogeography models.
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
animaColExt(Ext='crs', Col="dcr")
}

\Keyword{simulation}
\Keyword{Island Biogeography}

