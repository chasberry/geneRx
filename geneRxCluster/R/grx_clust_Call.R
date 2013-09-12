grx_clust_Call <-
function (chromoSts, chromoEnds, strt, grp, kvals, cutptExprs, 
    cutptFunExprs, tmpEnv, nperm) 
.Call("gRxC_cluster", chromoSts, chromoEnds, strt, grp, kvals, 
    cutptExprs, cutptFunExprs, tmpEnv, nperm, PACKAGE = "geneRxCluster")
