
##' Plot gRxCluster object clumps
##'
##' Plot Relative Frequencies of the two classes according to
##' region. Regions typically alternate between clusters and
##' non-clusters on each chromosome.
##' @title gRxPlotClumps
##' @param object result of gRxCluster
##' @param data (optional) GRanges like that from which args to gRxCluster were derived
##' @param seqlens (optional) seqlengths(data) or similar. Can be given if data is missing
##' @param panelExpr - an expression to evaluate after drawing each panel
##' @return NULL
##' @export
##' @example inst/ex-gRxPlotClumps.R
##' @include gRxCluster.R
##' @author Charles Berry
gRxPlotClumps <-
    function(object,data,seqlens,panelExpr=quote(grid()))
    {        
        if (missing(data)){
            dcall <- metadata(object)$call
            dcall[[1]] <- as.name("list")
            dcall <- dcall[1:4]
            data <- try(eval.parent(dcall))
            if (class(data)=="try-error")
                stop("Could not reconstruct data - Provide an explicit data argument")
            
            data <- 
                GRanges(seqnames=data$object,
                        IRanges(start=data$starts,
                                width=rep(1,length(data$starts))),
                        strand=rep("*",length(data$starts)),
                        group=data$group)
        }
        if (missing(seqlens)){
            seqlens <- seqlengths(object)
            if (any(is.na(seqlens))) seqlens <- try(seqlengths(data))
            if (class(seqlens)=="try-error" || any(is.na(seqlens))){
                message("Using highest starts as seglens")
                seqlens <- start(data)[cumsum(runLength(seqnames(data)))]
                names(seqlens) <- as.character(runValue(seqnames(data)))
            }
        }
        seqlengths(object)[names(seqlens)] <- seqlens
        
        max.len <- max(seqlens)
        p.null <- prop.table(table(data$group))[1]
        
        obj.gaps <- gaps(object)
        obj.gaps <- obj.gaps[strand(obj.gaps)=="*",]
        over.data <- findOverlaps(obj.gaps,data,ignore.strand=TRUE)            
        gap.tab <-
            table(factor(queryHits(over.data)),
                  data$group[subjectHits(over.data)])        
        
        mcols(obj.gaps)[as.numeric(rownames(gap.tab)),"value1"] <- gap.tab[,1]
        mcols(obj.gaps)[as.numeric(rownames(gap.tab)),"value2"] <- gap.tab[,2]
        
        gap.y <- obj.gaps$value1 / ( obj.gaps$value1 + obj.gaps$value2 )
        ungap.y <- object$value1 / ( object$value1 + object$value2 )
        
        npanels <- length(seqlens)
        par(mfrow = c( ceiling(npanels/3), 3 ),mar=c(0,0.5,0.2,0))
        for (i in names(seqlens)){
            gap.subset <- as.vector(seqnames(obj.gaps)==i)
            ungap.subset <- as.vector(seqnames(object)==i)
            xvals <- c(start(obj.gaps)[ gap.subset ],
                       start(object)[ ungap.subset ],
                       seqlens[i])
            xv.order <- order(xvals)
            yvals <- c( gap.y[gap.subset], 
                       ungap.y[ungap.subset],p.null)
            plot(xvals[xv.order],yvals[xv.order],ylim=0:1,xlim=c(-max.len/100,max.len),
                 xaxs="r",yaxs="i",
                 type='s',axes=F)
            box()
            segments( 1, p.null, seqlens[i], lty=3)
            text( max.len, 0.02, i, adj=c(1,0) )
            if (!is.null(panelExpr)) eval(panelExpr)
        }
    }
