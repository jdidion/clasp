#' Copy number analysis using the genoCN package.
#
#' genoCN was designed for copy number analysis of human tumor genotypes. Therefore, it
#' assumes that the default state is copy number 2 and all three genotypes possible (AA/AB/BB).
#' This is state 1 in both the genoCNV and genoCNA HMM. For inbred strains, the default state 
#' should actually be one of homozygosity (state 2).
#'
#' Short segments can be avoided by setting seg.nSNP to a larger number (default=3).
#'
#' The PFB format originated with PennCNV and is also used in genoCNA. A file is derived for each array.
#' We provide PFBs for MUGA and MegaMUGA. A future version of the software will include a function to
#' derive your own PFB, but for the time being please contact the author for help adpating this software
#' to a new array. A PFB has three columns: Chr, Position, PFB, with rownames being marker IDs.
#'
#' @param samples data.frame of sample info. At a minimum, name and type are required.
#' @param metrics list of data.frames (one per sample), each with copy number metrics (baf and lrr) 
#' for each marker.
#' @param pfb data.frame with data from PennCNV PFB file (see notes). 
#' @param ... additional arguments passed to .run.genoCNA.
#'
#' @return A list of genoCNAResult objects,
#'
#' @seealso \code{\link{estimate.copy.number}}
copy.number.genoCN <- function(samples, metrics, pfb, ...) {
    N <- nrow(samples)
    stopifnot(N==length(metrics))
    
    result <- lapply(1:N, function(i) {
        print(i)
        details <- .run.genoCNA(samples[i,"id"], metrics[[i]], pfb, inbred=samples[i,"type"] == "inbred", ...)
        
        summary <- do.call(rbind, by(details$ivls, details$ivls$chr, function(x) 
            tapply(x$end - x$start + 1, x$cn, sum)[as.character(0:4)]))
        summary[is.na(summary)] <- 0
        colnames(summary) <- 0:4
        summary <- summary / apply(summary, 1, sum)
        
        mean.cn <- as.vector(by(details$ivls, as.integer(details$ivls$chr), function(x) {
            sizes <- x$end - x$start + 1
            sum(sizes * x$cn) / sum(sizes)
        }))
        names(mean.cn) <- unique(details$ivls$chr)
        
        x <- list(mean.cn=mean.cn, summary=summary, details=details)
        class(x) <- "genoCNAResult"
        x
    })
    names(result) <- samples$id
    invisible(result)
}

.run.genoCNA <- function(name, metrics, pfb, params=NULL, args=NULL, inbred=FALSE, 
        wd=".", cleanup=FALSE, default.cn=c(2,2), cnv.only=FALSE) {
    # only use common markers
    i <- intersect(rownames(metrics), rownames(pfb))
    metrics <- metrics[i,]
    pfb <- pfb[i,]
    # remove markers with non-numeric BAF or LRR
    w <- .is.num(metrics$baf) & .is.num(metrics$lrr)
    metrics <- metrics[w,]
    pfb <- pfb[w,]
    # sort
    o <- order(pfb$Chr, pfb$Position)
    metrics <- metrics[o,]
    pfb <- pfb[o,]
    
    ## estimate initial transition probabilities
    
    # Load the default parameters
    data(init.Para.CNA, package="genoCN", envir=environment())

    # Modify the paramters for inbred samples
    if (inbred) {
        init.Para.CNA$trans.begin[1,] <- init.Para.CNA$trans.begin[1, c(2,1,3:9)]
    }
    # Define the transition probability matrix
    init.Para.CNA$trans.begin <- matrix(init.Para.CNA$trans.begin[1,], length(unique(pfb$Chr)), 9, byrow=TRUE)
    # Add any additional arguments
    if (!is.null(params)) {
        init.Para.CNA <- .merge.args(init.Para.CNA, params)
    }
    # Set the output dir
    curwd <- getwd()
    rmdir <- FALSE
    tryCatch({
        if (!exists(wd)) {
            created <- dir.create(wd, showWarnings=F, recursive=T)
            rmdir <- cleanup & created
        }
        setwd(wd)
        # run the algorithm
        if (cnv.only) {
            cnv.only <- rep(TRUE, nrow(metrics))
        }
        else {
            cnv.only <- NULL
        }
        genoCNA.args <- list(rownames(pfb), pfb$Chr, pfb$Position, metrics$lrr, metrics$baf, pfb$PFB, name, init.Para.CNA, 
            outputSeg=TRUE, outputSNP=0, contamination=TRUE, traceIt=0, cnv.only=cnv.only)
        if (!is.null(args)) {
            genoCNA.args <- .merge.args(genoCNA.args, args)
        }
        result <- do.call(genoCN::genoCNA, genoCNA.args)
        result$name <- name
        outfiles <- paste(name, c("segment.txt", "segment_viterbi.txt"), sep="_")
        result$segments <- read.table(outfiles[1], sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, comment.char='')
        if (is.null(default.cn)) {
            result$ivls <- result$segments[,1:5]
        }
        else {
            ivls <- NULL
            for (chrm in unique(pfb$Chr)) {
                chrm.ends <- range(pfb$Position[pfb$Chr==chrm])
                w <- result$segments$chr==chrm
                if (any(w)) {
                    iv <- result$segments[w, c(1:5)]
                    chrm.ivls <- iv
                    n <- nrow(iv)
                    if (iv[1,2] > chrm.ends[1]) {
                        chrm.ivls <- rbind(chrm.ivls, data.frame(
                            chr=chrm, start=chrm.ends[1], end=iv[1,2]-1, state=default.cn[1], cn=default.cn[2]))
                    }
                    if (iv[n,3] < chrm.ends[2]) {
                        chrm.ivls <- rbind(chrm.ivls, data.frame(
                            chr=chrm, start=iv[n,3]+1, end=chrm.ends[2], state=default.cn[1], cn=default.cn[2]))
                    }
                    if (n > 1) {
                        for (i in 1:(n-1)) {
                            if (iv[i,3] < iv[i+1,2]-1) {
                                chrm.ivls <- rbind(chrm.ivls, data.frame(
                                    chr=chrm, start=iv[i,3]+1, end=iv[i+1,2]-1, state=default.cn[1], cn=default.cn[2]))
                            }
                        }
                    }
                } 
                else {
                    chrm.ivls <- data.frame(chr=chrm, start=chrm.ends[1], end=chrm.ends[2], state=default.cn[1], cn=default.cn[2])
                }
                if (nrow(chrm.ivls) > 1) {
                    chrm.ivls <- .compress.intervals(chrm.ivls[order(chrm.ivls$start),], val.cols=c(4,5))
                }
                ivls <- rbind(ivls, chrm.ivls)
            }
            result$ivls <- ivls
        }
        if (cleanup) {
            for (f in outfiles) unlink(f)
            if (rmdir) unlink(wd)
        }
        result
    }, finally=setwd(curwd))
}

#' Plot the BAF, LRR and copy number states of a sample.
#'
#' @param metrics copy number metrics
#' @param result from calling copy.number.genoCN
#' @param cn.range the range of copy numbers to display
#' @param ... additional arguments to pass to \code{plot.metrics}.
copy.number.genoCN.plot <- function(sample.id, marker.set.id, metrics, 
        result, norm.params, db, cn.range=NULL, ...) {
    state.colors <- RColorBrewer::brewer.pal(9, "Set1")
    plot.cn.segments <- function(calls, metrics, other.data, chromosomes, xlim, ...) {
        segs <- other.data$ivls
        segs <- segs[as.character(segs[,1]) %in% names(chromosomes),]
        cn.range <- other.data$cn.range
        if (is.null(cn.range)) {
            cn.range <- range(segs[,5])
        }
        ylim <- c(cn.range[1] - 0.25, cn.range[2] + 0.25)
        plot(0:1, type="n", xlim=xlim, ylim=ylim, ylab="Copy Number", xaxt="n", yaxt="n", ...)
        axis(2, seq(cn.range[1],cn.range[2]), ...)
        segments(chromosomes[-1], ylim[1], chromosomes[-1], ylim[2], lty=2)
        start <- chromosomes[as.character(segs[,1])] + segs[,2]
        end <- chromosomes[as.character(segs[,1])] + segs[,3]
        alphas <- as.integer(segs[,1]) %% 2 == 0
        alphas[alphas==0] <- 0.5
        rect(start, segs[,5]-0.25, end, segs[,5]+0.25, col=rgb.alpha(state.colors[segs[,4]], alphas), border=NA)
    }
    metrics.plot(result$name, sample.id, marker.set.id, metrics, 
        norm.params, db, other.plots=list(plot.cn.segments), 
        other.data=list(list(ivls=result$details$ivls, cn=cn.range)), 
        ...)
}

# experimental
.est.purity <- function(result) {
    Beta.A.AB = result$mu.b[4,2]
    Beta.B.AB = result$mu.b[4,3]
     
    # the BAF for genotype AB
    Beta.AB = result$mu.b[1,2]
     
    # Update BAF estimates to take into account of systematic dye bias
    # this step is not necessary if there is no systematic bias
    Beta.A.AB = 0.5*Beta.A.AB/Beta.AB
    Beta.B.AB = 0.5 + 0.5*(Beta.B.AB - Beta.AB)/(1 - Beta.AB)
     
    # Re-estimate Beta.A.AB by averaging Beta.A.AB and 1 - Beta.B.AB
    Beta1 = 0.5*(Beta.A.AB + 1 - Beta.B.AB)
     
    # estimate tumor purity, since we examine the genotype A, nA = 1 and nB = 0
    nA = 1
    nB = 0
    pT = (1 - 2*Beta1)/(Beta1*(nA + nB - 2) + 1 - nB)
}
