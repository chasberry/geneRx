\name{cutpt.threshhold.binom}
\alias{cutpt.threshhold.binom}
\title{setup.threshhold.binom}
\usage{
  cutpt.threshhold.binom(k, n, val)
}
\arguments{
  \item{k}{the window widths}

  \item{n}{e.g. table(group)}

  \item{val}{maximum one tailed nominal alpha to determine
  cutpoints}
}
\value{
  a list
}
\description{
  find cutpoints
}
\details{
  find cutpoints based on expected number of False
  Discoveries in each tail. From gRxCluster, this is called
  after cutpt.filter.expr, so \code{x} - the sdiff object
  will be available.
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster}} for how and why this function is
  used
}

