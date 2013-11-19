\name{gRxPlot}
\alias{gRxPlot}
\title{gRxPlot}
\usage{
  gRxPlot(object, pi.0 = NULL,
    method = c("odds", "criticalRegions"), xlim = NULL,
    main = NULL, xlab = "log odds ratio",
    breaks = "Sturges", kvals = NULL, ...)
}
\arguments{
  \item{object}{either the results of
  \code{\link{gRxCluster}} or a list containing cutpoints
  for critical regions.}

  \item{pi.0}{the background proportion for vector 2}

  \item{method}{character vector of \dQuote{odds} and/or
  \dQuote{criticalRegions}}

  \item{xlim}{limits of the log odds histogram}

  \item{main}{a title for the panel(s)}

  \item{xlab}{label fgor the x-axis of the log odds plot}

  \item{breaks}{see \code{\link{hist}}}

  \item{kvals}{values to use in selecting a subset of the
  critical regions to display}

  \item{...}{other args to pass to the plotting routine(s)}
}
\value{
  see \code{\link{hist}}
}
\description{
  Plot Clumps and/or Critical Regions
}
\details{
  The results of a call to \code{\link{gRxCluster}} are
  plotted. Optionally, with \code{method="criticalRegions"}
  only the critical regions are plotted or with
  \code{method="odds"} the log odds only are plotted.
}
\examples{
x.seqnames <- rep(letters[1:3],each=500)
x.starts <- c(seq(1,length=500),seq(1,by=2,length=500),seq(1,by=3,length=500))
x.lens <- rep(c(5,10,15,20,25),each=20)
x.group <- rep(rep(c(TRUE,FALSE),length=length(x.lens)),x.lens)
## add a bit of fuzz:
x.group <- 1==rbinom(length(x.group),1,pr=ifelse(x.group,.8,.2))
x.kvals <- as.integer(sort(unique(x.lens)))
x.res <- gRxCluster(x.seqnames,x.starts,x.group,x.kvals)
gRxPlot(x.res)
rm( x.seqnames, x.starts, x.lens, x.group, x.kvals, x.res)
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxPlotClumps}} for a more fine grained
  display
}

