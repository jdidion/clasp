clasp
=====

Cell Line Authentication by SNP Profiling (CLASP).

This package supports the publication:

Didion JP et al, 'SNP array profiling of mouse cell lines identifies their strains 
    of origin and reveals cross-contamination and widespread 
    aneuploidy,' BMC Genomics 2014, 15:847, doi:10.1186/1471-2164-15-847.

# Installation

Installation requires the `devtools` package. Make sure not to build the vignette (set `build_vignettes=FALSE`, as shown below) unless you want to be waiting for a long time!

```r
install.packages(‘devtools’)
library(devtools)
devtools::install_github("jdidion/clasp/src/clasp", build_vignettes=FALSE)
```

# Usage

The vignette (vignettes/mouse-cell-lines.pdf) demonstrates how to use CLASP to reproduce the analyses from the paper. Because the vignette takes a very long time to run, it is not a part of the R package.
