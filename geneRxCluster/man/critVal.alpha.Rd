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
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster}} for how and why this function is
  used
}

