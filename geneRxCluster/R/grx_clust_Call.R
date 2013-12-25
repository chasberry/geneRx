grx_clust_Call <-
function (chromoSts, chromoEnds, strt, grp, kvals, cutptExprs, 
    cutptFunExprs, tmpEnv, nperm, sample_id, sample_tab) 
.Call("gRxC_cluster", chromoSts, chromoEnds, strt, grp, kvals, 
    cutptExprs, cutptFunExprs, tmpEnv, nperm, sample_id, sample_tab, 
    PACKAGE = "geneRxCluster")
