\name{threshhold.binom}
\alias{threshhold.binom}
\title{setup threshhold.binom}
\usage{
  threshhold.binom(k, val, p0)
}
\arguments{
  \item{k}{window widths}

  \item{val}{nominal alpha (length(val) %in%
  c(1,length(k)))}

  \item{p0}{null value}
}
\value{
  a list
}
\description{
  find cutpoints utility
}
\details{
  find cutpoints based on expected number of False
  Discoveries in each tail. Not for end-user use.
}
\author{
  Charles Berry
}

