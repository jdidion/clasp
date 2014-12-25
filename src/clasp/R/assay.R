#' Develop a new assay.
#'
#' Develop a new assay by selecting the most informative subset of \code{markers} for the
#' specified samples.
#'
#' An assay may be developed from a single marker set or from multiple intersecting marker sets. In
#' the latter case (i.e. if \code{markers} is a vector of two or more IDs), the intersection of all
#' the marker sets will be found using \code{intersect.marker.sets}.
#'
#' For the purposes of the consistency check and removing replicate samples, two samples are
#' considered replicates if they meet any of the following criteria:
#' 1. Identical names
#' 2. Identical non-NULL cell line names AND identical sources
#' 3. Null cell line names AND identical non-NULL backgrounds
#'
#' An SDP (sample distribution pattern) is the partitioning of a set of samples based on genotypes
#' for a single marker. It can be thought of as the phylogenetic tree with N braches, where N is the
#' number of possible genotypes. Consecutive markers with different SDPs are considered informative
#' even if their degree of linkage exceeds the specified \code{R^2} threshold.
#'
#' It is expected that all samples in the database have passed QC checks; this function does not
#' perform any sample-level QC.
#'
#' @param name Assay name
#' @param marker.set A vector of marker set IDs (see details) or a marker data.frame.
#' @param db The database connection.
#' @param sample.set A vector of sample set IDs, a sample data.frame or NULL if all samples should be used.
#' @param insert Whether or not to insert the new assay into the database.
#' @param description If insert=TRUE, the text description of the new marker set.
#' @param ob.name If insert=TRUE, the name of the new outbred-specific marker set.
#' @param ob.description If insert=TRUE, the text description of the new outbred-specific marker set.
#' @param exclude.cell.lines If TRUE, cell line samples are not used to develop the assay.
#' @param min.call.rate Minimum fraction of samples that have to have a non-missing (i.e. non-zero)
#' genotype (inclusive).
#' @param min.maf Minimum fraction of samples with non-missing genotype information that must have
#' a minor allele. The value is not inclusive, so a value of 0 means the MAF must be > 0, whereas a
#' value of NULL means the check is not performed.
#' @param max.hwe.pvalue The maximum p-value (inclusive) for a test of deviation from the expected
#' Hardy-Weinberg allele frequency.
#' @param consistency.check Whether a consistency check should be performed for replicates.
#' @param max.ld.r2 Maximum \code{R^2} value (inclusive) for pairwise linkage between consecutive markers.
#' @param max.ld.pvalue Maximum p-value (inclusive) for a test of pairwise linkage between
#' consecutive markers.
#' @param keep.informative.linked.markers If a linkage test is performed, this parameter specifies
#' that markers that exceed the \code{R^2} threshold should not be filtered if their SDPs differ
#' (see Details).
#' @param max.sdp.diff Maximum fraction by which two SDPs can differ and still be considered equal
#' for the purposes of \code{keep.informative.linked.markers}.
#' @param diag.coeffs genotype-speific coefficients for computing diagnostic values.
#' @param reduce Logical indicating whether the maximal marker set should be pruned using the
#' the following parameters.
#' @param min.diff When reduce=TRUE, the minimum allowable distance between any two markers.
#' @param min.markers When reduce=TRUE, the minimum number of markers in the final marker set. The
#' actual size of the final marker set may be smaller than this if the other filtering parameters
#' remove too many markers.
#' @param diag.threshold.iters A vector of thresholds to iterate through when pruining markers.
#'
#' @return A list with the following elements:
#' set.id: The ID of the new marker set created with the markers for this assay, or NULL if insert=FALSE.
#' markers: Full information for the set of markers included in the assay.
#' samples: Full information for the set of all samples considered for this assay.
#' sample.ixs: Indices of samples that were actually used to develop the assay (i.e. after removing
#' replicates).
#' replicates: A list  of replicate groups.
#' geno: MxN matrix of genotypes, where M is the assay markers and N is the assay samples.
#' diagnostic.values: MxN matrix of the diagnostic value for each SNP and for each sample.
#' distmat: NxN matrix of the pairwise number of unequal markers.
develop.assay <- function(name, marker.set, db, sample.set=NULL, insert=TRUE,
        description=NULL, ob.name=NULL, ob.description=NULL, exclude.cell.lines=TRUE,
		min.call.rate=0.8, min.maf=0, max.hwe.pvalue=NULL, consistency.check=TRUE, max.ld.r2=NULL,
		max.ld.pvalue=NULL, keep.informative.linked.markers=TRUE, max.sdp.diff=0.05, 
		diag.coeffs=list(N=1,AA=1,BB=1,AB=1), reduce=FALSE, min.diff=NULL, min.markers=NULL, 
		diag.threshold.iters=seq(0.25, 0.70, 0.05)) {
			
	marker.lookup.ids <- NULL
	sample.lookup.ids <- NULL
	
    # load marker and genotype information
	message("Loading markers...")
    if (length(marker.set) == 1) {
        markers <- get.markers(marker.set, db)
		marker.lookup.ids <- marker.set[1]
    }
    else {
        if (is.vector(marker.set)) {
            markers <- intersect.marker.sets(marker.set, db)
        }
        else {
            markers <- marker.set
        }
    }
    markers <- markers[order(markers$chromosome, markers$position_bp),]
	if (is.null(marker.lookup.ids)) {
		marker.lookup.ids <- markers$id
	}
	
    # load sample information
	message("Loading samples...")
    if (is.null(sample.set)) {
        samples <- get.all.samples(db)
    }
    else if (is.vector(sample.set)) {
		if (length(sample.set) == 1) {
			samples <- get.samples(sample.set, db, include.pedigree=T)
			sample.lookup.ids <- sample.set
		}
		else {
			samples <- do.call(rbind, lapply(sample.set, get.samples, db, include.pedigree=T))
		}
    }
    else {
        samples <- sample.set
    }
	if (exclude.cell.lines) {
		samples <- samples[is.na(samples$cell_line),]
	}
	if (is.null(sample.lookup.ids)) {
		sample.lookup.ids <- samples$id
	}
	
	# Load genotypes and order the rows and columns
	# the same as the marker and sample data frames
	message("Loading genotypes...")
    geno <- get.genotypes(marker.lookup.ids, db, sample.lookup.ids)
	ri <- match(markers$id, as.integer(rownames(geno)))
	ci <- match(samples$id, as.integer(colnames(geno)))
	geno <- geno[ri, ci]

    # Identify outbreds and treat them specially
	
    ob <- samples$type == "outbred"
    if (any(ob)) {
        ob.samples <- samples[ob,]
        samples <- samples[!ob,]
        ob.geno <- geno[,ob]
        geno <- geno[,!ob]
    }

    # Identify replicates. Replicates will be used to test markers for genotype consistency (if
    # consistency.check is non-null). Then, all but one from each set of replicates is removed
    # from the genotype matrix for the purpose of identifying diagnostic markers.

	message("Identifying replicate samples...")
    dup <- lapply(unique(samples$name[duplicated(samples$name)]), function(n) which(samples$name == n))
    names(dup) <- lapply(dup, function(x) as.character(min(x)))

    w.bg <- is.na(samples$cell_line) & !is.na(samples$background)
    if (any(w.bg)) {
        bg <- samples[w.bg, "background"]
        dup <- .merge.args(dup, sapply(unique(bg[duplicated(bg)]),
            function(b) which(w.bg & samples$background == b), simplify=FALSE))
    }

	if (!exclude.cell.lines) {
	    w.cl <- !is.na(samples$cell_line) & !is.na(samples$source)
	    if (any(w.cl)) {
	        lines <- samples[w.cl, c("cell_line","source")]
	        dup <- .merge.args(dup, apply(unique(lines[duplicated(lines),]), 1,
	            function(l) which(w.cl & samples$cell_line == l[1] & samples$source == l[2])))
	    }
	}
    
    has.ped <- which(!is.na(samples$mother_id) & !is.na(samples$father_id))

    # apply consistency check and compute error rate
	if (consistency.check) {
		message("Testing replicates for genotype consistency...")
        dup.consistent <- NULL
        for (i in 1:length(dup)) {
            before <- nrow(geno)
            geno <- geno[.consistent(geno[,dup[[i]]], FALSE, FALSE),]
            after <- nrow(geno)
            dup.consistent <- rbind(dup.consistent, c(i, before, after))
        }
        ped.consistent <- NULL
        for (i in has.ped) {
            before <- nrow(geno)
            geno <- geno[.f1.consistent(
                geno[,as.character(samples[i,c("id","mother_id","father_id")])], FALSE),]
            after <- nrow(geno)
            ped.consistent <- rbind(ped.consistent, c(i, before, after))
        }
    }

    # Save independent F1 genotypes
	f1.data <- NULL
	if (length(has.ped) > 0) {
		f1.dup.bg <- intersect(unique(samples[has.ped,"background"]), names(dup))
		f1.dup.idx <- which(names(dup) %in% f1.dup.bg)
		f1.keep <- sort(c(setdiff(has.ped, which(samples$background %in% f1.dup.bg)), 
			unlist(lapply(dup[f1.dup.idx], sample, 1))))
		f1.data <- list(samples=samples[f1.keep,], geno=geno[,f1.keep])
		dup <- dup[-f1.dup.idx]
	}
	
    # remove replicates and samples with pedigrees
	rem <- sort(unique(c(has.ped, unlist(dup))))
	rem <- setdiff(rem, unlist(lapply(dup, function(x) {
		# Make sure to retain a sample we want to impute, if any
		w <- samples[x,"impute"]
		if (any(w)) {
			x <- x[w]
		}
		if (length(x) > 1) {
			x <- sample(x, 1)
		}
		x
    })))
	
    geno <- geno[,-rem]
    samples <- samples[-rem,]

    af <- lapply(1:nrow(geno), function(i) .allele.freq(geno[i,]))

    # apply call rate filter
    if (!is.null(min.call.rate)) {
		message("Filtering by call rate...")
        w <- unlist(lapply(af, function(x) (1 - x$N.freq) >= min.call.rate))
        geno <- geno[w,]
        af <- af[w]
    }

    # apply MAF filter
    if (!is.null(min.maf)) {
		message("Filtering by MAF...")
        w <- unlist(lapply(af, function(x) min(x$p, x$q) > min.maf))
		geno <- geno[w,]
        af <- af[w]
    }

    # apply HWE filter
    if (!is.null(max.hwe.pvalue)) {
		message("Filtering by HWE...")
        w <- unlist(lapply(af, function(x) .hwe.test(x)$pvalue <= max.hwe.pvalue))
		geno <- geno[w,]
        af <- af[w]
    }

	markers <- markers[match(as.integer(rownames(geno)), markers$id),]

    # compute diagnostic value for each allele
	message("Computing diagnostic values...")
    allele.diag <- do.call(rbind, lapply(af, function(x) sapply(c("N","AA","BB","AB"),
        function(a) (1 - ((x[[a]] - 1) / x$n.geno)) * diag.coeffs[[a]])))
    # compute diagnostic value for each SNP
    snp.diag <- apply(allele.diag, 1, prod)

    # apply r2 filter
    # TODO: paralellize
    if (!is.null(max.ld.pvalue) || !is.null(max.ld.r2)) {
		message("Filtering by pairwise LD...")
		
	    marker.pair.ld <- do.call(rbind, lapply(unique(markers$chromosome), function(chrm) {
	        w.chrm <- which(markers$chromosome == chrm)
	        pairs <- data.frame(left=w.chrm[1:(length(w.chrm)-1)], right=w.chrm[2:length(w.chrm)])
	        cbind(pairs, do.call(rbind, apply(pairs, 1, function(i) {
	            ld <- .LD(af[unlist(i)])$ld
	            data.frame(r2=ld$r2, pvalue=ld$pvalue, linked=(
					!is.null(max.ld.pvalue) && ld$pvalue[1] < max.ld.pvalue) || (
					!is.null(max.ld.r2) && ld$r2 > max.ld.r2))
	        })))
	    }))
		
        if (any(marker.pair.ld$linked) && keep.informative.linked.markers) {
            # Test each pair to see if the SDP difference is above the threshold
            for (i in which(marker.pair.ld$linked)) {
                if (.sdp.diff(geno[unlist(marker.pair.ld[i,1:2]),]) > max.sdp.diff) {
                    marker.pair.ld[i,'linked'] <- FALSE
                }
            }
        }

        if (any(marker.pair.ld$linked)) {
            # identify intervals of consecutive linked markers
            n <- nrow(marker.pair.ld)
            w <- c(0, which(marker.pair.ld$linked[2:n] != marker.pair.ld$linked[1:(n-1)]), n)
            start <- w[1:(length(w)-1)] + 1
            end <- w[2:length(w)]

            rem <- NULL
            for (i in which(marker.pair.ld$linked[start])) {
                # if only two consecutive markers are linked, exclude the least informative one,
                # or pick one at random if they have the same diagnostic value
                if (start[i] == end[i]) {
                    ixs <- unlist(marker.pair.ld[start[i],1:2])
                    d <- snp.diag[ixs]
                    rem <- c(rem, ixs[sample(which(d == min(d)), 1)]) # TODO
                }
                # otherwise remove only the multi-linked markers
                else {
                    rem <- c(rem, unlist(marker.pair.ld[start[i],2]):unlist(marker.pair.ld[end[i],1]))
                }
            }
			
			markers <- markers[-rem,]
            geno <- geno[-rem,]
			af <- af[-rem]
			allele.diag <- allele.diag[-rem,]
			snp.diag <- snp.diag[-rem]
        }
    }

    # compute pairwise distances
    # TODO: paralellize
	message("Computing pairwise sample distances...")
	dists <- .pairwise.distances(geno)
	
    # if there are pairs with no differences at this point, there are either identical samples,
    # the filtering parameters are too stringent, or the marker set is not informative enough to
    # differentiate the current set of samples
    if (is.null(min.diff) || min.diff < 1) {
        min.diff <- 1
    }
    if (min(dists[,3]) < min.diff) {
        n <- samples$name
        stop(paste("The following pairs are have fewer than", min.diff, "differences:", paste(apply(
        	dists[which(dists[,3] < min.diff),], 1, function(x) paste(n[x[1]], n[x[2]], sep=" x ")), collapse="\n"), sep="\n"))
    }

    # select most informative markers
    if (reduce) {
		message("Reducing marker set...")
		
        # Calculate the theoretical minimum number of markers
		nsamples <- nrow(samples)
        min.markers <- max(min.markers, .estimate.min.markers(nsamples, 4, min.diff))

        nsnps <- nrow(geno)
        total.rem <- 0

        # This is the maximum diagnostic value that we can remove. We will slowly increment
        # this value when we can no longer remove any markers.
        for (diag.threshold in diag.threshold.iters) {
            if (nsnps <= min.markers) break

            threshold.rem <- 0

            for (pair in order(dists[,3])) {
                if (nsnps <= min.markers) break

                # For the current pair (starting with the most similar), determine which SNPs have
                # the lowest diagnostic value
                pair.diag <- apply(cbind(
                    .get.diagnostic.values(geno[,dists[pair, 1]], allele.diag),
                    .get.diagnostic.values(geno[,dists[pair, 2]], allele.diag)), 1, mean)

                while (TRUE) {
                    # number of SNPs to remove
                    below.threshold <- sum(pair.diag <= diag.threshold)
                    if (below.threshold == 0) {
                        break
                    }
                    # max we can remove before we hit the marker limit
                    max.rem1 <- nsnps - min.markers
                    if (max.rem1 <= 0) {
                        break
                    }
                    # max we can remove before we hit the pairwise distance limit
                    max.rem2 <- min(dists[,3]) - min.diff
                    # whether we have to check if this removal puts us below the threshold
                    check <- FALSE
                    if (max.rem2 == 0) {
                        num.rem <- 1
                        check <- TRUE
                    }
                    else {
                        num.rem <- min(below.threshold, max.rem1, max.rem2)
                    }

                    ixs <- order(pair.diag)[1:num.rem]
                    dist.diff <- sapply(1:nrow(dists),
                        function(i) sum(geno[ixs,dists[i,1]] != geno[ixs,dists[i,2]]))
                    temp.dist <- dists[,3] - dist.diff

                    if (check && min(temp.dist) < min.diff) {
                        # if we can't remove a SNP with going below one of the
                        # thresholds, move to the next pair
                        break
                    }

                    dists[,3] <- temp.dist
                    geno <- geno[-ixs,]
					markers <- markers[-ixs,]
                    allele.diag <- allele.diag[-ixs,]
                    pair.diag <- pair.diag[-ixs]
                    nsnps <- nsnps - length(ixs)
                    threshold.rem <- threshold.rem + num.rem

                    if (nsnps <= min.markers) break
                }
            }

            message(paste("Removed", threshold.rem, "SNPs at threshold", diag.threshold))
            total.rem <- total.rem + threshold.rem
        }
    }

	if (!is.null(f1.data)) {
		f1.data$geno <- f1.data$geno[rownames(geno),]
	}
	
	message("Getting diagnosic values and pairwise distances for final marker set...")
    diagnostic.values <- apply(geno, 2, .get.diagnostic.values, allele.diag)
    # NxN matrix of pairwise distances
	distmat <- .distance.matrix(dists, samples$name)

	# recompute LD
	message("(Re)computing pairwise marker LD on final marker set...")
    marker.pair.ld <- do.call(rbind, lapply(unique(markers$chromosome), function(chrm) {
        w.chrm <- which(markers$chromosome == chrm)
        pairs <- data.frame(left=w.chrm[1:(length(w.chrm)-1)], right=w.chrm[2:length(w.chrm)])
        cbind(pairs, do.call(rbind, apply(pairs, 1, function(i) {
            ld <- .LD(af[unlist(i)])$ld
            data.frame(r2=ld$r2, pvalue=ld$pvalue)
        })))
    }))

    # import the assay into the database
    set.id <- NULL
    if (insert) {
		message("Inserting assay marker set in database...")
        assay.marker.set.id <- import.marker.set(name, markers$id, db, description)
		assay.sample.set.id <- import.sample.set(name, db, description, samples$id)
    }

    # handle outbred samples
    ob.data <- NULL
	ob.assay.marker.set.id <- NULL
	ob.assay.sample.set.id <- NULL
    if (any(ob)) {
		message("Preparing assay for outbred samples...")
		
        ob.geno <- ob.geno[as.character(markers$id),]
		ob.bg <- ob.samples$background
        ob.dup <- duplicated(ob.bg)
		
        if (any(ob.dup)) {
            ob.consistent <- NULL
            for (b in unique(ob.bg[ob.dup])) {
                w <- which(ob.bg == b)
                before <- nrow(ob.geno)
                ob.geno <- ob.geno[.consistent(ob.geno[,w], FALSE, FALSE),]
                after <- nrow(ob.geno)
                ob.consistent <- rbind(ob.consistent, c(b, before, after))
            }
			ob.keep <- setdiff(1:nrow(ob.samples), unlist(lapply(unique(ob.bg[ob.dup]), function(x) {
				k <- w <- which(ob.bg==x)
				i <- ob.samples[w,"impute"]
				if (any(i)) {
					k <- k[i]
				}
				if (length(k) > 1) {
					k <- sample(k, 1)
				}
				setdiff(w, k)
            })))
            ob.geno <- ob.geno[,ob.keep]
			ob.samples <- ob.samples[ob.keep,]
        }
		
		ob.marker.ids <- as.integer(rownames(ob.geno))
        m <- match(ob.marker.ids, as.integer(rownames(geno)))
        ob.distmat <- .orthoganal.distance.matrix(ob.geno, geno[m,], 
			list(ob.samples$name, samples$name))

        ob.set.id <- NULL
        if (insert) {
            ob.assay.marker.set.id <- import.marker.set(ob.name, ob.marker.ids, db, ob.description)
			ob.assay.sample.set.id <- import.sample.set(ob.name, db, ob.description, ob.samples$id)
        }
        ob.data <- list(markers=markers[m,], samples=ob.samples, replicates=ob.dup, geno=ob.geno, distmat=ob.distmat)
    }

	message("Done!")

	assay <- data.frame(name=name, marker_set_id=assay.marker.set.id, sample_set_id=assay.sample.set.id,
		ob_name=ob.name, ob_marker_set_id=ob.assay.marker.set.id, ob_sample_set_id=ob.assay.sample.set.id,
        stringsAsFactors=FALSE)
	if (insert) {
		assay <- insert.assay(assay, db)
	}
	
	assay <- as.list(assay[1,])
	assay$inbred <- list(markers=markers, samples=samples, replicates=dup, geno=geno, 
			diagnostic.values=diagnostic.values, ld=marker.pair.ld, distmat=distmat)
	assay$f1 <- f1.data
	assay$outbred <- ob.data
	class(assay) <- "CLASPAssay"
	invisible(assay)
}

# Construct SDPs for each marker in a matrix of genotypes and then calculate all the pairwise
# differences. SDP similarity is measured by the sum of the pairwise shared samples across all
# genotypes divided by the total number of samples.
# TODO: paralellize
.sdp.diff <- function(g) {
    n <- nrow(g)
    exp <- ncol(g)
    sdps <- lapply(1:n, function(i) tapply(1:exp, unlist(g[i,]), identity))
    eq <- apply(data.frame(combn(n, 2)), 2, function(pair) {
        s1 <- sdps[[pair[1]]]
        s2 <- sdps[[pair[2]]]
        sum(sapply(as.character(0:3), function(i) length(intersect(s1[[i]],s2[[i]]))))
    })
    (exp - eq) / exp
}

# Use the Gilbert-Varshamov bound to estimate the minimum number of markers (n) that
# are required to construct a set of size x (i.e. number of samples) of sequences
# using an alphabet of size q (i.e. number of genotype classes) with minimum Hamming
# distance d between any two strings in the set.
.estimate.min.markers <- function(x, q, d) {
    if (q < 2) stop("q must be >= 2")
    if (d < 1) stop("d must be >= 1")

    A <- 0
    n <- max(1, floor(log(x, q))) # I am fairly certain by cannot prove that this is the minimum bound on n
    while (TRUE) {
        # calculate the Gilbert-Varshamov bound
        A <- (q^n) / sum(sapply(0:(d-1), function(j) choose(n, j) * ((q - 1)^j)))
        if (A < x) {
            n <- n + 1
        }
        else {
            return(n)
        }
    }
}

.get.diagnostic.values <- function(geno, diag) {
    result <- rep(NA, length(geno))
    for (a in 0:3) {
        w <- geno == a
        if (any(w)) {
            result[w] <- diag[w, a+1]
        }
    }
    result
}

#' Generate a summary of a CLASPAssay object.
#'
#' @param object CLASPAssay object.
#' @param ... unused.
#'
#' @return summary.CLASPAssay object.
summary.CLASPAssay <- function(object, ...) {
	value <- list()
	
	n <- nrow(object$inbred$markers)
	value$name <- object$name
	value$marker.set.id <- object$marker_set_id
	value$num.markers <- n
	value$sample.size <- nrow(object$inbred$samples)
	value$alignment.score <- 1 - (object$inbred$distmat[upper.tri(object$inbred$distmat)] / n)
	
	if (!is.null(object$f1)) {
		value$sample.size.f1 <- nrow(object$f1$samples)
	}
	
	no <- nrow(object$outbred$markers)
	value$name.ob <- object$ob_name
	value$marker.set.id.ob <- object$ob_marker_set_id
	value$num.markers.ob <- no
	value$sample.size.ob <- nrow(object$outbred$samples)
	value$alignment.score.ob <- 1 - (object$outbred$distmat[upper.tri(object$outbred$distmat)] / no)
	
	if ("diagnostic.values" %in% names(object$inbred)) {
		value$diagnostic <- apply(object$inbred$diagnostic.values,1,mean)
		value$ld.r2 <- object$inbred$ld$r2
	}
		
	class(value) <- "summary.CLASPAssay"
	value
}

#' Generate a string representation of a summary.CLASPAssay object.
#'
#' @param x an object of class summary.CLASPAssay
#' @param ... unused.
print.summary.CLASPAssay <- function(x, ...) {
	cat("\nSummary of assay ", x$name, ":\n", sep="")
	cat(rep("-", 18+nchar(x$name)), "\n", sep="")
	cat("Markers: Inbred =", x$num.markers, "; Outbred =", x$num.markers.ob, "\n")
	cat("Effective sample size: Inbred =", x$sample.size, "; Outbred =", x$sample.size.ob, "\n")
	cat("Pairwise alignment score: Inbred =", round(mean(x$alignment.score), 3), "(", 
        paste(round(range(x$alignment.score), 3), collapse=" - "), ") ; Outbred =", 
        round(mean(x$alignment.score.ob), 3), "(", 
        paste(round(range(x$alignment.score.ob), 3), collapse=" - "), ")\n")
	if ("diagnostic" %in% names(x)) { 
		cat("Mean marker diagnostic value: Inbred =", round(mean(x$diagnostic), 3), "\n")
		cat("Mean LD (r^2):", round(mean(x$ld.r2), 3))
	}
}

#' Execute an assay on a set of cell lines.
#' 
#' @param assay CLASPAssay object
#' @param db database connection
#' @param error.rate The error rate to use for PIA calculations.
#' @param cell.lines sample set IDs or vector of sample IDs of cell lines to test.
#' If null, all cell lines in the database will be tested.
#' @param impute.bgs boolean or vector of background (e.g. inbred strain) names. If
#' TRUE, synthetic F1s will be imputed for pairwise combinations of all backgrounds
#' with their impute field set to TRUE.
#' Only samples with impute column set to TRUE are used. If there is more than one
#' sample of a given background, one is selected randomly.
#' @param results.per.cl Number of matches to report for each cell line.
#' @param secondary.threshold match score below which a secondary background is searched.
#'
#' @return 
#' The return value is a list with the following elements:
#' 1. samples: Data frame of cell line sample data
#' 2. pairwise.matrix: Pairwise similarity matrix between all pairs of cell lines
#' 3. table: a summary data frame of the results with the following columns:
#'   id: Cell line sample ID
#'   primary.match: Name of best matching sample
#'   primary.match.score: Alignment score of primary match
#'   secondary.match: Name of best secondary match (if primary.match.score < secondary.threshold)
#'   secondary.match.score: Alignment score of secondary match
#'   pia: Probabilty of Incorrect Assignment
#' 4. details: a list the same length as the number of cell lines with the following elements:
#'   primary.matches
#'   secondary.matches
execute.assay <- function(assay, db, error.rate, cell.lines=NULL, impute.bgs=T, 
        results.per.cl=10, secondary.threshold=0.96) {
	stopifnot(class(assay) == "CLASPAssay")
	
	# ib = inbred background/control samples; ob = outbred background/control samples; cl = cell line samples
	
    message("Loading cell line annotations and genotypes...")
	if (is.null(cell.lines)) {
		cl.samples <- get.all.samples(db, where="cell_line IS NOT NULL")
		cell.lines <- cl.samples$id
	}
	else if (is.vector(cell.lines)) {
		cl.samples <- get.samples(cell.lines, db)
	}
	else {
		cl.samples <- cell.lines
		cell.lines <- cl.samples$id
	}
	
	cl.geno <- get.genotypes(assay$marker_set_id, db, cell.lines)
	
    # pairwise analysis of cell lines
    message("Computing pairwise identity of cell lines...")
    cl.pairs <- combn(1:nrow(cl.samples), 2)
    cl.eq <- apply(cl.pairs, 2, function(x) {
        sum(cl.geno[,x[1]] == cl.geno[,x[2]]) / nrow(cl.geno)
    })
    cl.sim <- .distance.matrix(cbind(t(cl.pairs), cl.eq), cl.samples$name)
	
	# compare cell lines to inbred controls
	ib.samples <- assay$inbred$samples
	ib.geno <- assay$inbred$geno
	if (!is.null(assay$f1)) {
		ib.samples <- rbind(ib.samples, assay$f1$samples)
		ib.geno <- cbind(ib.geno, assay$f1$geno)
	}

	controls <- ib.samples[,c("id","background","name")]
	all.names <- ib.samples$background
	
    if (impute.bgs != FALSE && any(ib.samples$impute)) {
        message("Imputing synthetic F1 samples...")
		w <- ib.samples$impute
		if (!is.logical(impute.bgs)) {
			w <- w & ib.samples$background %in% impute.bgs
		}
		im.ib.samples <- ib.samples[w,]
		im.ib.geno <- .force.matrix(ib.geno, 2, w)
		
        # Create synthetic F1s
        ib.imputed <- .impute.f1.genotypes(im.ib.geno)
        im.ib.pairs <- ib.imputed$sample.pairs
        im.ib.names <- paste("SYNTH_",
            controls[match(as.integer(im.ib.pairs[,1]), controls$id), 'name'], 'x',
            controls[match(as.integer(im.ib.pairs[,2]), controls$id), 'name'], sep='')
        colnames(ib.imputed$genotypes) <- im.ib.names
        all.names <- c(all.names, im.ib.names)
        ib.geno <- cbind(ib.geno, ib.imputed$genotypes)
    }
	
    # get the equality matrices for all cell lines
    # TODO: parallelize
    message("Comparing cell lines to control samples...")
    eq <- lapply(1:ncol(cl.geno), function(i) {
        cl.g <- cl.geno[,i]
        apply(ib.geno, 2, function(x) x == cl.g)
    })
    eq.frac <- lapply(eq, function(x) apply(x, 2, sum) / nrow(x))

	# compare cell lines to outbred controls, if any
    if (!is.null(assay$outbred)) {
        ob.geno <- assay$outbred$geno
        ob.samples <- assay$outbred$samples
        controls <- rbind(controls, ob.samples[,c('id','background','name')])
        all.names <- c(all.names, ob.samples$background)

        if (impute.bgs != FALSE && any(ob.samples$impute)) {
            message("Imputing synthetic F1s between outbred samples...")
			w <- ob.samples$impute
			if (!is.logical(impute.bgs)) {
				w <- w & ob.samples$background %in% impute.bgs
			}
			im.ob.samples <- ob.samples[w,]
			im.ob.geno <- .force.matrix(ob.geno, 2, w)
			
			# Create synthetic F1s from all combinations of inbreds AND outbreds,
			# only at markers selected for the outbred assay
			
            im.ob.ids <- as.character(im.ob.samples$id)
            im.ob.pairs <- do.call(rbind, lapply(im.ob.ids,
                function(x) cbind(rep(x, nrow(im.ib.samples)), im.ib.samples$id)))
            if (length(im.ob.ids) > 1) {
                im.ob.pairs <- rbind(im.ob.pairs, t(combn(im.ob.ids, 2)))
            }
            ob.imputed <- .impute.f1.genotypes(cbind(im.ib.geno[rownames(ob.geno),], im.ob.geno), im.ob.pairs)
            im.ob.names <- paste("SYNTH_",
                controls[match(as.integer(im.ob.pairs[,1]), controls$id), 'name'], 'x',
                controls[match(as.integer(im.ob.pairs[,2]), controls$id), 'name'], sep='')
            colnames(ob.imputed$genotypes) <- im.ob.names
            all.names <- c(all.names, im.ob.names)
            ob.geno <- cbind(ob.geno, ob.imputed$genotypes)
        }
        
        message("Comparing cell lines to outbred samples...")
        ob.eq <- lapply(1:ncol(cl.geno), function(i) {
            cl.g <- cl.geno[rownames(ob.geno), i]
            apply(ob.geno, 2, function(x) x == cl.g)
        })

        for (i in 1:length(eq.frac)) {
            eq.frac[[i]] <- c(eq.frac[[i]], apply(ob.eq[[i]], 2, sum) / nrow(ob.eq[[i]]))
        }
    }
    
    message("Computing primary match scores")
    matches <- lapply(eq.frac, function(x) {
        o <- order(x, decreasing=TRUE)[1:results.per.cl]
        data.frame(idx=o, name=all.names[o], score=x[o], stringsAsFactors=FALSE)
    })
    names(matches) <- cl.samples$name
	primary <- do.call(rbind, lapply(matches, function(x) x[1,]))
    
	# Secondary backgrounds
    find.secondary <- which(primary$score < secondary.threshold)
    secondary <- NULL
    if (length(find.secondary) > 0) {
        message(paste("Computing secondary match scores for", length(find.secondary), "cell lines"))
        secondary <- do.call(rbind, lapply(find.secondary, function(i) {
            c1 <- cl.geno[rownames(ob.geno), i]
            n2 <- matches[[i]][1,2]
            wn <- which(all.names == n2)
            if (wn <= ncol(ib.geno)) {
                c2 <- ib.geno[rownames(ob.geno),wn]
            }
            else {
                c2 <- ob.geno[,wn-ncol(ib.geno)]
            }
            ww <- c1==c2
            m1 <- apply(ib.geno[rownames(ob.geno),][!ww,], 2, function(x) sum(x == c1[!ww]))
            m2 <- apply(ob.geno[!ww,], 2, function(x) sum(x == c1[!ww]))
            m3 <- c(m1,m2)/sum(!ww)
            o <- order(m3, decreasing=TRUE)
            data.frame(cl=cl.samples[i,c(2,3,4,6,7)], matches[[i]][1,], secondary.idx=o[1], 
                secondary.name=all.names[o[1]], secondary.score=m3[o[1]], stringsAsFactors=FALSE)
        }))
        rownames(secondary) <- cl.samples[find.secondary,"name"]
    }
    
    # PIA
	message("Computing PIA values...")
    pia <- rep(NA, length(primary$name))
	for (m in unique(primary$name)) {
		i <- match(m, all.names)
		if (i <= ncol(ib.geno)) {
			gi <- ib.geno[,i]
			go <- ib.geno[rownames(ob.geno), i]
			hmax <- min(c(apply(ib.geno[,-i], 2, function(x) sum(gi != x)),
				  apply(ob.geno, 2, function(x) sum(go != x))))
		}
		else {
			i <- i - ncol(ib.geno)
			g <- ob.geno[,i]
			hmax <- min(c(apply(ib.geno[rownames(ob.geno),], 2, function(x) sum(g != x)),
  			      apply(ob.geno[,-i], 2, function(x) sum(g != x))))
		}
		pia[m == primary$name] <- hmax
	}
	pia <- Rmpfr::mpfr(error.rate, 100) ^ pia
    
    summary <- data.frame(
        id=cl.samples$id,
        primary.match=primary$name,
        primary.match.score=primary$score,
        secondary.match=NA,
        secondary.match.score=NA,
        pia=sprintf("%.1E",pia),
        row.names=cl.samples$name,
        stringsAsFactors=F)
    summary[find.secondary, 4:5] <- secondary[,10:11]
    
    details <- lapply(1:nrow(cl.samples), function(i) {
        n <- cl.samples$name[i]
        list(primary.matches=matches[[n]], secondary.matches=secondary[[n]], pia=pia[i])
    })
    
    result <- list(samples=cl.samples, pairwise.matrix=cl.sim, summary=summary, details=details)
    class(result) <- "CLASPAssayResults"
    
    message("Done!")
    
    invisible(result)
}

#' Compute copy-number metrics.
#'
#' Compute copy-number metrics (BAF and LRR) for cell line samples based on
#' intensity values.
#'
#' Computing LRR requires a large number of reference samples. We provide reference parameters 
#' for the MUGA and MegaMUGA array. A future version of the software will provide tools to 
#' derive new reference parameters. In the meantime, please contact the author for help adapting
#' new arrays to work with this software.
#'
#' @param marker.set.id ID of marker set containing markers that match params
#' @param params reference parameters (see note)
#' @param db Database connection.
#' @param cell.lines Sample set ID for cell line samples, or NULL to select all cell lines.
#' @param platform Platform name.
#' @param normalize whether to normalize the intensity values.
#'
#' @return list of metrics (one per sample)
compute.metrics <- function(marker.set.id, params, db, cell.lines=NULL, platform="Illumina", normalize=TRUE) {
    if (is.null(cell.lines)) {
    	cl.samples <- get.all.samples(db, where="cell_line IS NOT NULL")
    }
    else {
    	cl.samples <- get.samples(cell.lines, db)
    }
    cl.data <- get.platform.data(marker.set.id, platform, db, cell.lines, as.matrix=TRUE, format=2)
    .compute.metrics <- get(paste("compute.metrics", platform, sep="."))
    cl.metrics <- .compute.metrics(cl.data, params, normalize)
    invisible(cl.metrics)
}

#' Create a DilutionSeries object from a set of samples in the database.
#'
#' A DilutionSeries is used to predict sample contamination. To create a
#' dilution series, make in vitro mixtures of two geneticially divergent
#' samples in a series of ratios (the more the better). Then genotype these
#' samples and import them into the database.
#'
#' For MUGA and MegaMUGA arrays, you can use the pre-computed DilutionSeries
#' \code{default.dilution.series}.
#'
#' Ratios table columns:
#' 1. sample.id
#' 2. control The proportion of the mixture containing the control sample
#' 3. contaminant The proportion of the mixture containing the contaminant sample
#' control.frac and contam.frac must be numeric values summing to 1.0 for each row.
#'
#' @param marker.set.id ID of marker set for which to compute dilution series metrics.
#' @param sample.set.id ID of sample set containing dilution series samples.
#' @param ratios A data frame with three columns (see note).
#' @param params Reference parameters.
#' @param db Database connection
#' @param platform array platform.
#'
#' @return a DilutionSeries object.
compute.dilution.series <- function(marker.set.id, sample.set.id, ratios, params, 
        db, platform="Illumina") {
    samples <- get.samples(sample.set.id, db, include.pedigree=FALSE)
    data <- get.platform.data(marker.set.id, platform, db, sample.set.id)
    .compute.metrics <- get(paste("compute.metrics", platform, sep="."))
    series <- list(N=nrow(params), ratios=ratios, metrics=.compute.metrics(data, params))
    class(series) <- "DilutionSeries"
    invisible(series)
}

#' Estimate contamination using platform-specific (i.e. intensity) information.
#'
#' @param cl.metrics Metrics for cell lines.
#' @param dilution.series A DilutionSeries, or a path to an RData file containing a
#' pre-computed DilutionSeries.
#' @param baf.thresholds vector of two BAF thresholds (homozygous and heterozygous).
#'
#' Estimating BAF thresholds requires a large number of heterozygous control samples.
#' For our paper, we used the F1 samples in the database to compute the mean mirrored
#' BAF for homozygous and heterozygous calls and rounded up and down, respectively, to 
#' two decimal places. A future version of this software will provide a more robust
#' method of determining the thresholds, as well as tools for generating reference
#' parameters for new arrays. Until then, please contact the author for help adapting
#' new arrays to work with this software.
#'
#' @return
#' Numeric vector with names as sample IDs and values as predictions as the fraction
#' of contaminant cells.
estimate.contamination <- function(cl.metrics, dilution.series, baf.thresholds=c(0.02, 0.46)) {
    if (is.character(dilution.series)) {
        dilution.series <- load.objects(dilution.series)
    }
    x <- dilution.series$ratios$contaminant
    y <- do.call(rbind, lapply(dilution.series$metrics, .shift, baf.thresholds))
    n <- y[,1] / y[,2]
    ss <- smooth.spline(x, y[,3] * n)
    ypred <- predict(ss, seq(0,1,0.01))
    # only use the increasing part of the line for prediction
    w <- max(which(ypred$y == (max(ypred$y))))
    ypredx <- ypred$x[1:w]
    ypredy <- ypred$y[1:w]
    
    cl.y <- do.call(rbind, lapply(cl.metrics, .shift, baf.thresholds))
    cl.n <- unlist(lapply(1:nrow(cl.y), function(i) cl.y[i,1] / cl.y[i,2]))
    cl.x <- unlist(lapply(cl.y[,3] * cl.n, function(y) {
        d <- abs(y - ypredy)
        ypredx[which(d==min(d))]
    }))
    names(cl.x) <- names(cl.metrics)
    return(list(contamination=cl.x,
        ypred=ypred, ds=data.frame(x=x, y=y[,3] * n), cl=data.frame(x=cl.x, y=cl.y[,3] * cl.n)))
}

.shift <- function(x, baf.thresholds) {
    b <- x$baf
    n <- !is.na(b)
    z <- pmin(b[n], 1-b[n])
    w <- z <= 0.25
    z[w] <- z[w] - baf.thresholds[1]
    z[!w] <- baf.thresholds[2] - z[!w]
    data.frame(N=sum(z>0), D=sum(n), z=mean(z[z > 0]))
}

#' Estimate copy number for a set of cell lines.
#'
#' @param samples data.frame of sample info. At a minimum, name and type are required.
#' @param metrics Metrics for the cell line samples.
#' @param cn.args List of rguments for the copy number estimation function.
#' @param method Copy number estimation method.
#'
#' @return list of copy number predictions. Each list item is a list
#' with three elements:
#' 1. mean.cn Vector of numeric values, the weighted mean copy number for
#' each chromosome.
#' 2. summary Table with one or more entries for each chromosome. Each
#' entry gives the fraction of the chromosome that has a given copy number.
#' 3. details Original results from the specific copy number prediction
#' method that was used.
estimate.copy.number <- function(samples, metrics, cn.args, method="genoCN") {
    func <- get(paste("copy.number", method, sep="."))
    invisible(do.call(func, c(list(samples, metrics), cn.args)))
}