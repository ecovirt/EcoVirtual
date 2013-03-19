\Name{extGame}

\Title{
Zero-sum game
}

\Description{
Simulates a zero-sum game between two competitors in a finite number of individuals.
}

\Usage{
extGame(bet, total, tmax)
}

\Arguments{
  \item{bet}{
    bet size of each competitor on each time.
}
  \item{total}{
    total amount of 'utility' in the game.
}
} \item{tmax}{
    maximum game time (number of cicles).
}
}

\Details{
In game theory and economic theory, a zero–sum game is a mathematical representation of a situation in which a participant's gain (or loss) of utility is exactly balanced by the losses (or gains) of the utility of the other participant(s). If the total gains of the participants are added up, and the total losses are subtracted, they will sum to zero. A zero–sum game is also called a strictly competitive game while non-zero–sum games can be either competitive or non-competitive.
}

\Value{
'extGame' returns a graph with the number of individuals of each competitor on each cicle.

'extGame' also returns an invisible vector with the results of the loser on each bet.
}

\References{

}

}
\Author{
Alexandre Adalardo de Oliveira <aleadalardo@gmail.com>
}

\See Also{
simHub, randWalk
}

\Examples{
extGame(bet=1,total=20)
extGame(bet=1,total=100)

}

\Keyword{simulation}
\Keyword{neutral theory}
