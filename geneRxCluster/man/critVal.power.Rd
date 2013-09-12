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
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster}} for how and why this function is
  used
}

