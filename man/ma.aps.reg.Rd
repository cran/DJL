\name{ma.aps.reg}
\alias{ma.aps.reg}
\title{Combinatorial search (all possible subset) for regression analysis}
\description{Implements combinatorial (exhaustive) search algorithm, aka all-possible-subsets regression. As opposed to the sequential approach (stepwise, forward addition, or backward elimination) that has a potential bias resulting from considering only one variable for selection at a time, all possible combinations of the independent variables are examined, and sets satisfying designated conditions are returned.}
\usage{ma.aps.reg(dv, iv, min=1, max, mad=FALSE, aic=FALSE, bic=FALSE,
           model.sig=TRUE, coeff.sig=TRUE, coeff.vif=TRUE, coeff.cor=FALSE)}
\arguments{
  \item{dv}{Dependent variable (\emph{r} by \emph{1} matrix)}
  \item{iv}{Independent variable(s) (\emph{r} by \emph{c} matrix)}
  \item{min}{Minimum number of independent variable to explore (>=\emph{1})}
  \item{max}{Maximum number of independent variable to explore (<=\emph{r}/10)}
  \item{mad}{Returns mean absolute deviation when \code{TRUE}}
  \item{aic}{Returns Akaike's information criterion when \code{TRUE}}
  \item{bic}{Returns Bayesian information criterion when \code{TRUE}}
  \item{model.sig}{Returns models statistically significant only when \code{TRUE}}
  \item{coeff.sig}{Returns models with statistically significant coefficients only when \code{TRUE}}
  \item{coeff.vif}{Returns models with allowable level of multicollinearity only when \code{TRUE}}
  \item{coeff.cor}{Returns models without suppression effects only when \code{TRUE}}
}
\references{Hair, Joseph F., et al. Multivariate data analysis. Vol. 7. \emph{Upper Saddle River}, NJ: Pearson Prentice Hall, 2006.}
\author{Dong-Joon Lim, PhD}
\examples{
# Load airplane dataset
data(dataset.airplane.2017)

# ready
dv<-dataset.airplane.2017[,2,drop=FALSE]
iv<-dataset.airplane.2017[,3:7,drop=FALSE]
  
# go
ma.aps.reg(dv,iv,1,3,mad=TRUE,coeff.cor=TRUE)
}