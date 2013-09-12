
##' best contiguous region
##'
##' prune each end of the region using loglik criterion
##' @title prune.loglik
##' @param x a GRanges object
##' @param p.null the probability of category 1 (FALSE)
##' @return genomewide probability of the FALSE category
##' @author Charles Berry
prune.loglik <- function(x,p.null=0.5)
{
    x.best <- function(x){
        abs.depth <- abs(x$depth)
        high.scores <- range(which(abs.depth==max(abs.depth)))
        GRanges(
            seqnames = seqnames( x )[ 1 ] ,
            IRanges(start=start(x)[high.scores[1]],
                    end=end(x)[high.scores[2]]),
            depth=x$depth[high.scores[1]],
            clump.id=x$clump.id[high.scores[1]])
    }
    
    llfun <- function(x,y)
        ifelse(x>0,x*log(x),0)+ifelse(y>0,y*log(y),0)-(x+y)*log(x+y)
    llnull <- function(x,y)
        x*log(p.null)+y*log(1-p.null)
    
    
    x.max <- function(res.elt,x.elt){
        tmp1 <- subjectHits(findOverlaps(ranges(res.elt),ranges(x.elt)))
        ##    tmp1 <- subjectHits(findOverlaps((res.elt),(x.elt)))
        tab <- as.data.frame(mcols(x.elt[,2:3]))
        ctab <- colSums( tab[ unique(tmp1),] )
        
        ltab <- rbind( tab[(1:nrow(tab))<min(tmp1),],
                      ctab)
        
        rtab <- rbind(ctab,
                      tab[(1:nrow(tab))>max(tmp1),])
        
        rtab[["c1"]] <- cumsum(rtab[,1])
        rtab[["c2"]] <- cumsum(rtab[,2])
        rtab[["rest1"]] <- rev(cumsum(rev(rtab[,1])))-rtab[,1]
        rtab[["rest2"]] <- rev(cumsum(rev(rtab[,2])))-rtab[,2]
        r.llk <- llnull(rtab$rest1, rtab$rest2 ) + llfun(rtab$c1, rtab$c2 )
        ltab[["c1"]] <- rev(cumsum(rev(ltab[,1])))
        ltab[["c2"]] <- rev(cumsum(rev(ltab[,2])))
        ltab[["rest1"]] <- cumsum(ltab[,1])-ltab[,1]
        ltab[["rest2"]] <- cumsum(ltab[,2])-ltab[,2]
        l.llk <- llnull(ltab$rest1, ltab$rest2 ) + llfun(ltab$c1, ltab$c2 )
        l.max.at <- min(which(l.llk==max(l.llk)))
        r.max.at <- max(tmp1) - 1 + max(which(r.llk==max(r.llk)))
        res <- GRanges(
            seqnames = seqnames( x.elt )[ 1 ] ,
            IRanges(start=start(x.elt)[l.max.at],
                    end=end(x.elt)[r.max.at]))
        mcols(res) <-
            t(as.matrix(c(colSums(tab[l.max.at:r.max.at,]),
                          clump.id=as.vector(x.elt$clump.id[1]))))
        
        res
    }
    x <- split(x,x$clump.id)
    res.list <- lapply(x,x.best)
    clmp.list <- mapply(x.max,res.list,as.list(x),SIMPLIFY=FALSE)
    unlist(GRangesList( clmp.list ))
}

##' critical region cutpoints
##'
##' This version uses alpha and will find TFD
##' @title critical regions
##' @param k - window width(s)
##' @param p0 - length 2 probabilities
##' @param alpha - two tailed 
##' @param posdiff - position difference matrix
##' @return list of cutoffs and attributes
##' @seealso \code{\link{gRxCluster}} for how and why this function is used
##' @author Charles Berry
##' @export
critVal.alpha <-
    function(k,p0,alpha,posdiff){
        ns <- colSums(!is.na(posdiff))
        tails <-
            mapply(function(kelt,n){
                mat <- cbind(
                    low=pbinom(0:kelt,kelt,p0[2]),
                    hi=pbinom(-1+0:kelt,kelt,p0[2],lower.tail=F)
                    )
                low_cut <- tail(which(mat[,'low']<alpha/2),1)
                high_cut <- head(which(mat[,'hi']<alpha/2),1)
                res <- 
                    c(low= if (length(low_cut)) low_cut else 0,
                    up= if (length(high_cut)) high_cut-1 else nrow(mat))
                mat <- cbind(2*n*mat,mat)
                colnames(mat) <- c("target.low","target.hi","low","hi")
                attr(res,"fdr") <- mat
                attr(res,"target") <- alpha*n
                res
            },k,ns,SIMPLIFY=FALSE)
        tails}

##' critical region cutpoints
##'
##' This version uses TFD and will find alpha implicitly
##' @title critical regions
##' @param k  window width(s)
##' @param p0  length 2 probabilities
##' @param target - two tailed 
##' @param posdiff - position difference matrix
##' @param ns the number of windows passing filter at each k
##' @return list of cutoffs and attributes
##' @seealso \code{\link{gRxCluster}} for how and why this function is used
##' @author Charles Berry
##' @example inst/ex-critVal-target.R
##' @export
critVal.target <-
    function(k,p0,target,posdiff=NULL,ns){
        p0 <- prop.table(p0) # just in case raw counts were passed
        stopifnot(length(p0)==2)
        if (!is.null(posdiff)) ns <- colSums(!is.na(posdiff))
        if (length(k) != length(ns))
            stop("length(k) must = ncol(posdiff) or length(ns)")
        tails <-
            mapply(function(kelt,n){
                mat <- cbind(
                    low=pbinom(0:kelt,kelt,p0[2]),
                    hi=pbinom(-1+0:kelt,kelt,p0[2],lower.tail=F)
                    )
                low_cut <- tail(which(mat[,'low']<target/n/2),1)
                high_cut <- head(which(mat[,'hi']<target/n/2),1)
                res <- 
                    c(low= if (length(low_cut)) low_cut else 0,
                    up= if (length(high_cut)) high_cut-1 else nrow(mat))
                mat <- cbind(2*n*mat,mat)
                colnames(mat) <- c("target.low","target.hi","low","hi")
                attr(res,"fdr") <- mat
                attr(res,"target") <- target 
                res
            },as.list(k),as.list(ns),SIMPLIFY=FALSE)
        tails}

##' critical region cutpoints
##'
##' This version uses power and TFD and will limit windows screened
##' @title critical regions
##' @param k - window width(s)
##' @param p0 - length 2 probabilities
##' @param target - false discoveries wanted
##' @param pwr - desired power
##' @param odds - alternative odds ratio
##'  @return list of cutoffs and attributes
##' @seealso \code{\link{gRxCluster}} for how and why this function is used
##' @author Charles Berry
##' @export
critVal.power <-
  function(k,p0,target,pwr=0.8,odds=7)
{     
    pi.0 <- prop.table(p0)[2]
    pi.alt <- plogis( qlogis( pi.0 ) + c(-1,1) * log(odds) )
    ctpts <-
        as.data.frame(
            t(matrix(
                sapply(1:2,
                       function(x) (x==1) +
                       qbinom(1-pwr,k,pi.alt[x],lower.tail=x==2)),
                ncol=2)))
    names(ctpts) <- k
    target.min <-
        lapply(k,function(x){
            y <- dbinom(0:x,x,pi.0)
            cbind(low=cumsum(y),hi=rev(cumsum(rev(y))))})    
    alpha <-
        mapply(function(x,y)
               max(abs(pbinom(x-1,y,pi.0)-0:1)),
               ctpts,
               k)
    target.nj <- target/alpha/2
    res <- mapply(function(x,y,z){
        attr(x,"fdr") <- cbind( target= 2 * y * trunc(z),y )
        colnames(attr(x,"fdr")) <- c("target.low","target.hi","low","hi")
        attr(x,"n") <- z
        attr(x,"target") <- target
        x},
                  ctpts,target.min,target.nj,
                  SIMPLIFY=FALSE)
    attr(res,"filter.fun") <-
        function(x) {
            n <- target.nj
            sapply(seq_along(n),
                   function(j) 
                   quantile(x[,j],
                            min(n[j]/sum(!is.na(x[,j])),1.0),
                           na.rm=TRUE))}
    res
            
}

##' find cutpoints
##'
##' find cutpoints based on expected number of False Discoveries in
##' each tail. From gRxCluster, this is called after
##' cutpt.filter.expr, so \code{x} - the sdiff object will be
##' available.
##' @title setup.threshhold.binom
##' @param k    the window widths   
##' @param n    e.g. table(group)                           
##' @param val  maximum one tailed nominal alpha to determine cutpoints
##' @return a list
##' @seealso \code{\link{gRxCluster}} for how and why this function is used
##' @author Charles Berry
cutpt.threshhold.binom <-
    function(k,n,val)
    {
        if (val<0.0 || val>1.0) stop("val must be in [0,1]") 
        p.nought <- prop.table(n)[1]
        threshhold.binom(k,val,p.nought)
    }

##' find cutpoints utility
##'
##' find cutpoints based on expected number of False Discoveries in
##' each tail. Not for end-user use.
##' @title setup threshhold.binom
##' @param k  window widths 
##' @param val  nominal alpha (length(val) %in% c(1,length(k)))
##' @param p0 null value
##' @return a list
##' @author Charles Berry
threshhold.binom <-
    function(k,val,p0)
    {
        ##  ccnt matrix of counts of efun_vals
        ##  k the window widths for each column of ccnt
        ##  n table(efun_vals)
        ##  val the target for alpha

        if (length(val)==1) val <- rep(val,length(k))
        stopifnot(length(val)==length(k))
        
        bino.prs <-
            lapply(k, function(ki) {
                dh <- dbinom(ki:0, ki, p0 )
                dim(dh) <- c(1,ki+1)
                colnames(dh) <- 0:ki
                dh})
        bino.median.index <-
            sapply(bino.prs,
                   function(x) which(cumsum(x)>=0.5)[1])
        e.tails <-
            lapply(bino.prs,
                   function(hp){
                       lower <- cumsum(hp)
                       upper <- rev(cumsum(rev(hp)))
                       cbind(up=upper,low=lower)})
        for (i in seq_along(e.tails))
            attr(e.tails[[i]],"val") <- val[i] 
        
        ctv.result <- function(et){
            
            fdr_low <- et[,'low']
            fdr_up <- et[,'up']
            low_cut <- tail(which(fdr_low < attr(et,"val")),1)
            high_cut <- head(which(fdr_up < attr(et,"val")),1)
            res <- 
                c(low= if (length(low_cut)) low_cut else 0,
                  up= if (length(high_cut)) high_cut-1 else length(fdr_up))
            attr(res,"fdr") <- et[,c("low","up")]
            attr(res,"val") <- attr(et,"val")
            res
        }
                
        lapply( e.tails,ctv.result)
    }

##' Plot a set of cutpoints
##'
##' not for users. not exported
##' @title plot.cutpoints
##' @param crit - a cutpoint object see \code{\link{gRxCluster}}
##' @param pi.0 - optional null value to plot
##' @param kvals - which cutpoints to includein the plot
##' @param ... passed to barplot
##' @return list with components of \dQuote{bar.x} (the value of
##' \code{hist()}), \dQuote{kvals} (window widths plotted), and
##' \dQuote{pi.0} (the input vlue of \code{pi.0}
##' @author Charles Berry
plot.cutpoints <-
    function(crit, pi.0=NULL, kvals=NULL, ...){
        tmp.k <- 
            sapply(crit,function(x) nrow(attr(x,"fdr"))-1)
        both.bars <- apply(rbind(0,sapply(crit,c),1+tmp.k),2,diff)
        if (!is.null(kvals)){
            k.indx <- match(kvals,tmp.k)
            stopifnot(all(is.finite(k.indx)))
            tmp.k <- tmp.k[k.indx]
            both.bars <- both.bars[, k.indx, drop=F]
        }

        bar.x <- barplot(both.bars,col=c("gray","white","gray"),
                         xlab=expression("Sites in Window" == w[j]),names.arg=tmp.k,
                         ...
                         )
        if (!is.null(pi.0))
            points(bar.x,pi.0*tmp.k,pch=19,cex=1.0)
        invisible(list(bar.x=bar.x,kvals=tmp.k,pi.0=pi.0))
    }
