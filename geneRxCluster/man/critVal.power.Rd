\name{critVal.power}
\alias{critVal.power}
\title{critical regions}
\usage{
  critVal.power(k, p0, target, pwr = 0.8, odds = 7)
}
\arguments{
  \item{k}{- window width(s)}

  \item{p0}{- length 2 probabilities}

  \item{target}{- false discoveries wanted}

  \item{pwr}{- desired power}

  \item{odds}{- alternative odds ratio}
}
\value{
  list of cutoffs and attributes
}
\description{
  critical region cutpoints
}
\details{
  This version uses power and TFD and will limit windows
  screened
}
\examples{
# symmetric odds:
crit <-
    critVal.power(5:25,c(1,1),5,pwr=0.8,odds=7)
crit[[1]]
sapply(crit,c)
# 5:1 odds
asymmetric.crit <-
    critVal.power(5:25,c(1,5),5,pwr=0.8,odds=7)
# show the critical regions
par(mfrow=c(1,2))
gRxPlot(crit,method="critical")
gRxPlot(asymmetric.crit,method="critical")
rm(crit,asymmetric.crit)
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster}} for how and why this function is
  used
}

