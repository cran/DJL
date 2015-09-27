\name{target.spec.dea}
\alias{target.spec.dea}
\title{Spec target setting using DEA}
\description{Employs \code{inverse DEA} to estimate the specifications to achieve a determined efficiency.}
\usage{target.spec.dea(xdata, ydata, date, t, dt, dmu, et="c", alpha=NULL, beta=NULL,
                wv, rts, sg="ssm", ftype="d", ncv=NULL, env=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
  \item{t}{A vantage point from which the RoC is captured}
  \item{dt}{Delta \emph{t} i.e., specs are estimated within PPS at \emph{t+dt}}
  \item{dmu}{DMU whose inputs(or outputs) are to be estimated}
  \item{et}{Efficiency target; default value ("c") retains the current efficiency}
  \item{alpha}{Perturbed input(s) of designated DMU (\emph{1} by \emph{m} matrix)}
  \item{beta}{Perturbed output(s) of designated DMU (\emph{1} by \emph{s} matrix)}
  \item{wv}{Weight vector for scalarization (\emph{1} by \emph{m} or \emph{s} matrix)}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization} \cr
    \code{"min"  }{Date-sum maximization}
  }
  \item{ftype}{Frontier type \cr
    \code{"d"  }{Dynamic frontier (default)} \cr  
    \code{"s"  }{Static frontier}
  }  
  \item{ncv}{Non-controllable variable index for internal NDF (\emph{1} by \emph{(m+s)} binary matrix)}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1} matrix)}
}
\value{
  \item{\code{$alpha}}{Estimated input(s)}
  \item{\code{$beta}}{Estimated output(s)}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
}
\references{
Lim, Dong-Joon, "Inverse DEA with frontier changes for new product target setting." \emph{Submitted} (2015).
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{roc.dea}} RoC calculation using DEA \cr
  \code{\link{target.arrival.dea}} Arrival target setting using DEA \cr
  \code{\link{target.spec.dea}} Spec target setting using DEA
}
\examples{
# Reproduce Table 4 in Lim, D-J. (2015)
  # Load engine dataset
    data(dataset.engine.2015)
  
  # Subset for forced induction systems
    fis<-subset(dataset.engine.2015,grepl("^.C..",dataset.engine.2015[,8]))
  
  # ready
    # Suppose one wants to estimate Porsche 911 turbo s' engine specs
    # to retain its current competitiveness with downsized 3.5 litre engine in 2018.
    # What might be the minimum specs to achieve this goal
    # considering the technological changes we've seen so far?
    # Plus, the CEO wants to put more emphasis on the torque improvement than HP.
      x<-matrix(as.numeric(fis[,4]),ncol=1)
      y<-matrix(as.numeric(fis[,6:7]),ncol=2)
      d<-matrix(as.numeric(fis[,2]),ncol=1)
      a<-matrix(c(3.5),ncol=1)
      w<-matrix(c(0.3,0.7),ncol=2)
    
  # go
    target.spec.dea(x,y,d,2015,3,262,alpha=a,wv=w,rts="vrs",sg="min")
}