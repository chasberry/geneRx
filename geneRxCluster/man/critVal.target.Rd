\name{critVal.target}
\alias{critVal.target}
\title{critical regions}
\usage{
  critVal.target(k, p0, target, posdiff = NULL, ns)
}
\arguments{
  \item{k}{window width(s)}

  \item{p0}{length 2 probabilities}

  \item{target}{- two tailed}

  \item{posdiff}{- position difference matrix}

  \item{ns}{the number of windows passing filter at each k}
}
\value{
  list of cutoffs and attributes
}
\description{
  critical region cutpoints
}
\details{
  This version uses TFD and will find alpha implicitly
}
\examples{

# symmetric odds:
crit <- critVal.target(5:25,c(1,1),1,ns=rep(10,21))
crit[[1]]
sapply(crit,c)
# 5:1 odds
asymmetric.crit <- critVal.target(5:25,c(1,5),1,ns=rep(10,21))
# show the critical regions
par(mfrow=c(1,2))
gRxPlot(crit,method="critical")
gRxPlot(asymmetric.crit,method="critical")
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster}} for how and why this function is
  used
}

