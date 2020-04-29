\name{plp}
\alias{plp}
\title{Print LP object}
\description{Print an LP object line by line.}
\usage{plp(x)}
\arguments{
  \item{x}{LP object defined by make.lp function in lpSolve library}
}
\references{
Berkelaar, Michel, Kjell Eikland, and Peter Notebaert. "lpsolve: Open source (mixed-integer) linear programming system." \emph{Eindhoven U. of Technology} 63 (2004).
}
\author{Dong-Joon Lim, PhD}
\examples{
# Declare an LP object
lp.temp <- make.lp(0, 61)

# Print the LP
plp(lp.temp)
}