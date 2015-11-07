\name{roc.sf}
\alias{roc.sf}
\title{Rate of change (RoC) calculation using SF}
\description{Employs \code{dm.sf} over time to calculate RoCs. This function is valid only when multiplicative form of directional vector is used.}
\usage{roc.sf(xdata, ydata, date, t, rts, g, w=NULL, sg="ssm", ftype="d")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
  \item{t}{A vantage point from which the RoC is captured}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)} matrix)}
  \item{w}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s} matrix)} 
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization} \cr
    \code{"min"  }{Date-sum minimization}
  }
  \item{ftype}{Frontier type \cr
    \code{"d"  }{Dynamic frontier (default)} \cr  
    \code{"s"  }{Static frontier}
  }
}
\value{
  \item{\code{$eff_r}}{Efficiency at release (i.e., at each production date)}
  \item{\code{$eff_t}}{Efficiency at \code{t}}
  \item{\code{$lambda_t}}{Intensity vector at \code{t}}
  \item{\code{$eft_date}}{Effective date}
  \item{\code{$roc_past}}{RoC observed from the obsolete DMUs in the past}
  \item{\code{$roc_avg}}{Average RoC}
  \item{\code{$roc_local}}{Local RoC}
}
\references{
D.-J. Lim,  Internal combustion engine race: naturally aspirated vs turbo/super-charged, \emph{working paper} (2015).
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{roc.dea}} RoC calculation using DEA \cr
  \code{\link{target.arrival.dea}} Arrival target setting using DEA \cr
  \code{\link{target.spec.dea}} Spec target setting using DEA
}
\examples{
# Reproduce Mercedes-Benz CLA45 AMG's local RoC in Table 5 in Lim, D-J. (2015)
  # Load engine dataset
    data(dataset.engine.2015)
  
  # Parameters
    c<-matrix(as.numeric(dataset.engine.2015[,3]),ncol=1)
    x<-matrix(as.numeric(dataset.engine.2015[,4]),ncol=1)
    y<-matrix(as.numeric(dataset.engine.2015[,5:7]),ncol=3)
    d<-matrix(as.numeric(dataset.engine.2015[,2]),ncol=1)
    g<-matrix(c(rep(0,nrow(x)),y),nrow=nrow(x),ncol=4)
    w<-matrix(c(1,0,0),ncol=3)

  # Subset for 4 cylinder engines
    x_s<-subset(x,c==4)
    y_s<-subset(y,c==4)
    d_s<-subset(d,c==4)
    g_s<-subset(g,c==4)

  # Calc local Roc
    roc.sf(x_s,y_s,d_s,2014,"crs",g_s,w,"min")$roc_local[348,]
}