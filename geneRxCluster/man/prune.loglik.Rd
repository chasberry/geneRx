\name{prune.loglik}
\alias{prune.loglik}
\title{prune.loglik}
\usage{
  prune.loglik(x, p.null = 0.5)
}
\arguments{
  \item{x}{a GRanges object}

  \item{p.null}{the probability of category 1 (FALSE)}
}
\value{
  genomewide probability of the FALSE category
}
\description{
  best contiguous region
}
\details{
  prune each end of the region using loglik criterion
}
\author{
  Charles Berry
}

