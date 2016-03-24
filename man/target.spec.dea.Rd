\name{target.spec.dea}
\alias{target.spec.dea}
\title{Spec target setting using DEA}
\description{Employs \code{inverse DEA} to estimate specifications(in/out-puts) to achieve a predetermined efficiency.}
\usage{target.spec.dea(xdata, ydata, date=NULL, t=NULL, dt=NULL, dmu, et="c",
                alpha=NULL, beta=NULL, wv, rts, sg="ssm", ftype="d",
                ncv=NULL, env=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{t}{A vantage point from which the RoC is captured}
  \item{dt}{Delta \emph{t} i.e., specs are estimated within PPS at \emph{t+dt}}
  \item{dmu}{DMU whose inputs(or outputs) are to be estimated}
  \item{et}{Efficiency target; default value ("c") retains the current efficiency}
  \item{alpha}{Perturbed input(s) of designated DMU (\emph{1} by \emph{m})}
  \item{beta}{Perturbed output(s) of designated DMU (\emph{1} by \emph{s})}
  \item{wv}{Weight vector for scalarization (\emph{1} by \emph{m} or \emph{s})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization} \cr
    \code{"min"  }{Date-sum minimization}
  }
  \item{ftype}{Frontier type \cr
    \code{"d"  }{Dynamic frontier (default)} \cr  
    \code{"s"  }{Static frontier}
  }  
  \item{ncv}{Non-controllable variable index(binary) for internal NDF (\emph{1} by \emph{(m+s)})}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1})}
}
\value{
  \item{\code{$alpha}}{Estimated input(s)}
  \item{\code{$beta}}{Estimated output(s)}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
}
\references{
Lim, Dong-Joon, "Inverse DEA with frontier changes for new product target setting." 
\emph{European Journal of Operational Research} (In press) (2016).

Wei, Quanling, Jianzhong Zhang, and Xiangsun Zhang. "An inverse DEA model for inputs/outputs estimate." \emph{European Journal of Operational Research} 121.1 (2000): 151~163.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{roc.dea}} RoC calculation using DEA \cr
  \code{\link{target.arrival.dea}} Arrival target setting using DEA
}
\examples{
# Reproduce Example 2 in Wei, Q. et al.(2000)
  # ready
  x<-matrix(c(1,1,1),3)
  y<-matrix(c(4,8,5,8,4,5),3)
  a<-matrix(1.8,1)
  w<-matrix(c(0.5,0.5),1)
  
  # go
  target.spec.dea(x,y,dmu=3,alpha=a,wv=w,rts="crs")$beta

# Reproduce Table 4 in Lim, D-J. (2016)
  # Load engine dataset
    data(dataset.engine.2015)
  
  # Subset for forced induction systems
    fis<-subset(dataset.engine.2015,grepl("^.C..",dataset.engine.2015[,8]))
  
  # ready
    # Suppose one wants to estimate Porsche 911 turbo s' engine specs
    # to retain its current competitiveness with downsized 3.5 litre engine in 2018.
    # What might be the minimum specs to achieve this goal
    # considering the technological changes we've seen so far?
    # Plus, the CEO wants to put more emphasis on the torque improvement over HP.
      d<-subset(fis,select=2)
      x<-subset(fis,select=4)
      y<-subset(fis,select=6:7)
      a<-as.matrix(3.5)
      w<-matrix(c(0.3,0.7),1)
    
  # go
    target.spec.dea(x,y,d,2015,3,262,alpha=a,wv=w,rts="vrs",sg="min")$beta
}