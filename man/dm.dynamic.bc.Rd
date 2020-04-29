\name{dm.dynamic.bc}
\alias{dm.dynamic.bc}
\title{Dynamic DEA in the presence of intertemporal Budget Constraints}
\description{Employs the Farrell measure on carry-over budget as well as input or output}
\usage{dm.dynamic.bc(xdata, ydata, zdata, bdata, rts="crs", orientation="i", wv=NULL)}
\arguments{
  \item{xdata}{Input array (\emph{n} by \emph{m} by \emph{t})}
  \item{ydata}{Output array (\emph{n} by \emph{s} by \emph{t})}
  \item{zdata}{Budget(spent) array (\emph{n} by \emph{b} by \emph{t})}
  \item{bdata}{Budget(secured) array (\emph{n} by \emph{b})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{orientation}{Orientation of the measurement \cr
    \code{"i"  }{Input-orientation (default)} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{wv}{Weight vector for scalarization (\emph{1} by \emph{m} or \emph{s})}
}
\value{
  \item{\code{$eff.s}}{System Efficiency}
  \item{\code{$eff.t}}{Period Efficiency}
  \item{\code{$lambda}}{Intensity vectors}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
  \item{\code{$zslack}}{Budget(spent) slack}
  \item{\code{$aslack}}{Budget(available) slack}
}
\references{
Lim, D.-J., M.-S., Kim, & K.-W., Lee. (2020). "A revised dynamic DEA model with budget constraints." \emph{International Transactions in Operational Research (In press)}.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA
}
\examples{
# Load data
  df.io <- array(c(2, 4, 8, 4, 1, 2, 2, 2, 3, 6, 12, 6,
                   5, 4, 3, 8, 1, 1, 1, 1, 5, 4,  3, 8),
                 c(4, 3, 2), 
                 dimnames = list(LETTERS[1:4], c("X", "Y", "z"), c("t1", "t2")))
  df.Z.0 <- array(c(9, 12, 18, 24), c(4, 1), dimnames = list(LETTERS[1:4], c("Z^0")))

# Run
  dm.dynamic.bc(df.io[,1,], df.io[,2,], df.io[,3,], df.Z.0)
}