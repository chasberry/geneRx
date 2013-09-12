
##' Plot Clumps and/or Critical Regions
##'
##' The results of a call to \code{\link{gRxCluster}} are plotted. Optionally, with
##' \code{method="criticalRegions"} only the critical regions
##' are plotted or with \code{method="odds"} the log odds only are plotted.
##' @title gRxPlot
##' @param object either the results of \code{\link{gRxCluster}} or a
##' list containing cutpoints for critical regions.
##' @param pi.0 the background proportion for vector 2
##' @param method character vector of \dQuote{odds} and/or \dQuote{criticalRegions}
##' @param xlim limits of the log odds histogram  
##' @param main a title for the panel(s)
##' @param xlab label fgor the x-axis of the log odds plot
##' @param breaks see \code{\link{hist}}
##' @param kvals values to use in selecting a subset of the critical
##' regions to display
##' @param ... other args to pass to the plotting routine(s)
##' @seealso \code{\link{gRxPlotClumps}} for a more fine grained display
##' @return see \code{\link{hist}}
##' @include gRxCluster.R
##' @export
##' @author Charles Berry
gRxPlot <- function(object,pi.0=NULL, method=c("odds","criticalRegions"),
                    xlim=NULL,main=NULL,xlab="log odds ratio",
                    breaks="Sturges",kvals=NULL,...){
    match.methods <- match.arg(method, several.ok=TRUE)

    if (any(match.methods=="odds" )){
        if (!inherits(object,"DataFrame"))
            stop("need Granges object for method='odds'")
        vals <- as.data.frame(mcols(object)[,c("value1","value2")])
        log.odds <- qlogis((vals$value2+0.5)/(rowSums(vals)+1)) - 
            if (is.null(pi.0)) qlogis((sum(vals$value2)+0.5)/(sum(vals)+1))
            else
                qlogis( pi.0)
        hx <- hist( log.odds , breaks=breaks, plot=FALSE )
        louter <- if (is.null(xlim)) c(-1,1)*max(abs(hx$breaks)) else xlim
        plot(hx,xlim=louter,xlab=xlab,main=main,...)
    }

    if (any(match.methods=="criticalRegions")){
        if (inherits(object,"GRanges"))
            ctpt <- metadata(object)$criticalValues
        if (is.list(object)) ctpt <- object
        plot.cutpoints(ctpt,pi.0=pi.0,kvals=kvals,main=main, ...)    
    }
}
