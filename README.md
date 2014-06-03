geneRx
======

Tools for Data on Revtroviral Vectors as used in gene therapy  

# Installing

`geneRxCluster` is on Bioconductor. The easiest (and usually best)
method to install it is to issue these commands in your `R` session:

    source("http://bioconductor.org/biocLite.R")
    biocLite("geneRxCluster")

# Files in this repo

The files found here are possibly out of date. They follow this format:

-   **geneRxCluster:** development version of the geneRxCluster R package

-   **geneRxCluster_x.y-z.tar.gz:** version x.y-z of the package as tar.gz
    source - ready to install as source

-   **geneRxCluster_x.y-z.zip:** version x.y-z of the package as Windows
    binary installation

-  **geneRxCluster_x.y-z.tgz:**  version x.y-z of the package as Mac OS X
    binary installation

The package requires the GenomicRanges and IRanges bioConductor
packages. To install them, issue the commands in Installing and they
will automatically be installed. For manual installations, issues
these commands in R


    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicRanges")


