\name{gRxPlotClumps}
\alias{gRxPlotClumps}
\title{gRxPlotClumps}
\usage{
  gRxPlotClumps(object, data, seqlens,
    panelExpr = quote(grid()))
}
\arguments{
  \item{object}{result of gRxCluster}

  \item{data}{(optional) GRanges like that from which args
  to gRxCluster were derived}

  \item{seqlens}{(optional) seqlengths(data) or similar.
  Can be given if data is missing}

  \item{panelExpr}{- an expression to evaluate after
  drawing each panel}
}
\value{
  NULL
}
\description{
  Plot gRxCluster object clumps
}
\details{
  Plot Relative Frequencies of the two classes according to
  region. Regions typically alternate between clusters and
  non-clusters on each chromosome.
}
\examples{

x.seqnames <- rep(letters[1:3],each=50)
x.starts <- c(seq(1,length=50),seq(1,by=2,length=50),seq(1,by=3,length=50))
x.lens <- rep(c(5,10,15,20,25),each=2)
x.group <- rep(rep(c(TRUE,FALSE),length=length(x.lens)),x.lens)
## add a bit of fuzz:
x.group <- 1==rbinom(length(x.group),1,pr=ifelse(x.group,.8,.2))
x.kvals <- as.integer(sort(unique(x.lens)))
x.res <- gRxCluster(x.seqnames,x.starts,x.group,x.kvals)
gRxPlotClumps(x.res)
rm( x.seqnames, x.starts, x.lens, x.group, x.kvals, x.res)
}
\author{
  Charles Berry
}

