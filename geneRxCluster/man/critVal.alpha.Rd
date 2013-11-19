\name{critVal.alpha}
\alias{critVal.alpha}
\title{critical regions}
\usage{
  critVal.alpha(k, p0, alpha, posdiff)
}
\arguments{
  \item{k}{- window width(s)}

  \item{p0}{- length 2 probabilities}

  \item{alpha}{- two tailed}

  \item{posdiff}{- position difference matrix}
}
\value{
  list of cutoffs and attributes
}
\description{
  critical region cutpoints
}
\details{
  This version uses alpha and will find TFD
}
\examples{
# symmetric odds:
crit <- critVal.alpha(5:25,c(1,1)/2,alpha=0.05,
                      matrix(1,nr=50,nc=21))
crit[[1]]
sapply(crit,c)
# 5:1 odds
asymmetric.crit <- critVal.alpha(5:25,c(1,5)/6,alpha=0.05,
                      matrix(1,nr=50,nc=21))
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

