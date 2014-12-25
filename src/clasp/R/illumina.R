#' Compute metrics using Illumina array data.
#'
#' @param data list of matricies of platform-specific data (one per variable)
#' @param params reference parameters
#' @param normalize whether to normalize intensity data
#'
#' @return list of metrics (baf and lrr), one per sample.
compute.metrics.Illumina <- function(data, params, normalize=TRUE) {
    N <- ncol(data[[1]])
    metrics <- lapply(1:N, function(i) {
        intens <- data.frame(X=data$X[,i], Y=data$Y[,i], row.names=rownames(data$X))
        
        # HACK: The MUGA and MegaMUGA parameters were calculated based on data in which
        # the X and Y values were swapped for some markers. For that reason, there is 
        # an extra column in the paramters called "switched." We swap the x and y values
        # for each row in which switched=TRUE.
        if ("switched" %in% names(params) && any(params$switched)) {
            n <- names(params$switched)[params$switched]
            m <- match(n, rownames(intens))
            temp <- intens[m, "X"]
            intens[m, "X"] <- intens[m, "Y"]
            intens[m, "Y"] <- temp
        }
        
        if (normalize) {
            metrics <- .tQN(intens, params$clusters, params$QN.thresholds)
            metrics <- .wave.correction(metrics, params$gcmodel)
        }
        else {
            metrics <- .default.copy.number.metrics(intens, params$baflrr)
        }
    })
    names(metrics) <- colnames(data[[1]])
    invisible(metrics)
}

# This method produces non-normalized intensities and is
# mainly intended for comparison purposes.
.default.copy.number.metrics <- function(intens, params) {
    m <- match(rownames(params), rownames(intens))
    params <- params[!is.na(m),]
    intens <- intens[m[!is.na(m)],]

    valid <- intens$X >= 0 & intens$Y >= 0
    params <- params[valid,]
    intens <- intens[valid,]

    theta <- (2 / pi) * atan(intens$Y / intens$X)
    theta[is.nan(theta)] <- 0
    r <- intens$X + intens$Y

    baf <- rep(NaN, length(theta))
    r.exp <- rep(NaN, length(theta))
    s.exp <- rep(NaN, length(theta))

    # BAF is the linear interpolation of the sample theta value between the two nearest
    # canonical clusters. It is clipped to the range [0,1].
    het.missing <- is.na(params$thetaAB)
    if (any(het.missing)) {
        idx1 <- het.missing & theta <= params$thetaAA
        idx2 <- het.missing & theta > params$thetaAA & theta < params$thetaBB
        idx3 <- het.missing & theta >= params$thetaBB

        baf[idx1] <- 0
        baf[idx2] <- ((theta[idx2] - params$thetaAA[idx2]) / (params$thetaBB[idx2] - params$thetaAA[idx2]))
        baf[idx3] <- 1

        r.exp[idx1] <- params$rAA[idx1]
        r.exp[idx2] <- (params$mrA[idx2] * theta[idx2]) + params$brA[idx2]
        r.exp[idx3] <- params$rBB[idx3]

        s.exp[idx1] <- params$sAA[idx1]
        s.exp[idx2] <- (params$msA[idx2] * theta[idx2]) + params$bsA[idx2]
        s.exp[idx3] <- params$sBB[idx3]
    }

    idx1 <- !het.missing & theta <= params$thetaAA
    idx2 <- !het.missing & theta > params$thetaAA & theta < params$thetaAB
    idx3 <- !het.missing & theta == params$thetaAB
    idx4 <- !het.missing & theta > params$thetaAB & theta < params$thetaBB
    idx5 <- !het.missing & theta >= params$thetaBB

    baf[idx1] <- 0
    baf[idx2] <- 0.5 * ((theta[idx2] - params$thetaAA[idx2]) / (params$thetaAB[idx2] - params$thetaAA[idx2]))
    baf[idx3] <- 0.5
    baf[idx4] <- 0.5 + (0.5 * ((theta[idx4] - params$thetaAB[idx4]) / (params$thetaBB[idx4] - params$thetaAB[idx4])))
    baf[idx5] <- 1

    # LRR is the log transformation of the ratio between observed and expected intensity.
    # Expected intensity is calculated by an interpolation of theta along the line between
    # the centroids of the two nearest intenstype clusters.
    r.exp[idx1] <- params$rAA[idx1]
    r.exp[idx2] <- (params$mrA[idx2] * theta[idx2]) + params$brA[idx2]
    r.exp[idx3] <- params$rAB[idx3]
    r.exp[idx4] <- (params$mrB[idx4] * theta[idx4]) + params$brB[idx4]
    r.exp[idx5] <- params$rBB[idx5]

    lrr <- log2(r / r.exp)

    # EXPERIMENTAL: investigating the use of Z-scores rather than/in addition to LRR.
    s.exp[idx1] <- params$sAA[idx1]
    s.exp[idx2] <- (params$msA[idx2] * theta[idx2]) + params$bsA[idx2]
    s.exp[idx3] <- params$sAB[idx3]
    s.exp[idx4] <- (params$msB[idx4] * theta[idx4]) + params$bsB[idx4]
    s.exp[idx5] <- params$sBB[idx5]

    z <- (r - r.exp) / s.exp

    valid <- !(is.nan(baf) | is.nan(r.exp) | is.nan(z))
    data.frame(baf=baf[valid], lrr=lrr[valid], z=z[valid], row.names=rownames(intens)[valid], stringsAsFactors=FALSE)
}

# Perform tQN normalization of intensity data.
# Adapted from code provided by Johan Staaf.
# Staaf et al., BMC Bioinformatics, 2008, doi:10.1186/1471-2105-9-409.
.tQN <- function(intens, clusters, QN.thresholds) {
    # Remove any rows that aren't in the cluster file
    i <- intersect(rownames(intens), rownames(clusters))
    ref <- clusters[i,]
    intens <- intens[i,]
    
    # Quantile normalization
    inorm <- normalizeQuantiles(intens)
    inorm$X[inorm$X < 0] <- 0
    inorm$Y[inorm$Y < 0] <- 0
    
    # Calculate R
    R <- inorm$X + inorm$Y

    # Thesholding of QN effect
    aff.x <- which((inorm$X / intens$X) > QN.thresholds[1])
    inorm$X[aff.x] <- QN.thresholds[1] * intens$X[aff.x]

    aff.y <- which((inorm$Y / intens$Y) > QN.thresholds[2])
    inorm$Y[aff.y] <- QN.thresholds[2] * intens$Y[aff.y]

    T <- .theta(inorm)
    
    # Calculate tQN X and Y to fit theta and R
    Y <- R * tan(T * pi/2) / (1 + tan(T * pi/2))
    X <- R - Y

    # Calculate BAF and LRR from corrected R and T
    baflrr <- .baf.lrr(R, T, ref)

    data.frame(baflrr, X=X, Y=Y, row.names=rownames(intens))
}

.baf.lrr <- function(R, T, ref) {
    BAF <- T
    LRR <- R
    
    med.tAA <- median(ref$AA_T_Mean, na.rm=TRUE)
    med.tBB <- median(ref$BB_T_Mean, na.rm=TRUE)
    
    for (i in 1:length(T)) {
        if (!.is.num(T[i])) {
            BAF[i] <- NaN
            LRR[i] <- NaN
            next
        }
    
        th <- T[i]
        rr <- R[i]
        rAA <- ref$AA_R_Mean[i]
        rAB <- ref$AB_R_Mean[i]
        rBB <- ref$BB_R_Mean[i]
        tAA <- ref$AA_T_Mean[i]
        tAB <- ref$AB_T_Mean[i]
        tBB <- ref$BB_T_Mean[i]

        e.tAA <- .is.num(tAA)
        e.tAB <- .is.num(tAB)
        e.tBB <- .is.num(tBB)

        # 0: Test for inconsistencies between tAA/tAB/tBB and rAA/rAB/rBB
        if (((e.tAA & e.tAB) && tAA > tAB) ||
                ((e.tAA & e.tBB) && tAA > tBB) ||
                ((e.tAB & e.tBB) && tAB > tBB)) {
            BAF[i] <- LRR[i] <- NaN
        }
        # 1: Triple blank SNP
        else if (!(e.tAA|e.tAB|e.tBB)) {
            BAF[i] <- LRR[i] <- NaN
        }
        # 2: Blank for AB, AA, while positive for BB
        else if (!(e.tAA|e.tAB) & e.tBB) {
            if (th >= tBB) {
                BAF[i] <- 1
                LRR[i] <- ifelse(rBB <= 0, NaN, log2(rr / rBB))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 3: Blank for AB, BB, while positive for AA
        else if (e.tAA & !(e.tAB|e.tBB)) {
            if (th <= tAA) {
                BAF[i] <- 0
                LRR[i] <- ifelse(tAA <= 0, NaN, log2(rr / tAA))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 4: Blank for AB while positive for AA & BB
        else if (e.tAA & !e.tAB & e.tBB) {
            # no AB cluster exist for this SNP, while AA & BB exists. Set it to the closest of AA or BB
            min.index <- which.min(c(abs(tAA - th), abs(tBB - th)))
            if (min.index == 1 && th < tAA) {
                BAF[i] <- 0
                LRR[i] <- ifelse(tAA <= 0, NaN, log2(rr / tAA))
            }
            else if (min.index != 1 && th >= tBB) {
                BAF[i] <- 1
                LRR[i] <- ifelse(rBB <= 0, NaN, log2(rr / rBB))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 5: Blank for AA while positive for AB & BB
        else if (!e.tAA & e.tAB & e.tBB) {
            if (th >= tBB) {
                BAF[i] <- 1
            }
            # 5.1: SNP is "correctly between" ref$AB_T_Mean and ref$BB_T_Mean
            else if (th >= tAB) {
                # interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
                BAF[i] <- 0.5 + 0.5 * (th - tAB) / (tBB - tAB)  
            }
            # 5.2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be 
            # between ref$AA_T_Mean and ref$AB_T_Mean where it normally should not NOT BE.
            else {
                BAF[i] <- ifelse(th < med.tAA, 0, 0.5 * (th - med.tAA) / (tAB - med.tAA))             
            }
            
            if (th > tAB) {
                eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
                LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
            }
            else {
                LRR[i] <- NaN
            }
        }
        # 6: Blank for BB while positive for AA & AB
        else if (e.tAA & e.tAB & !e.tBB) {
            if (th < tAA) {
                BAF[i] <- 0
            }
            # 6.1: SNP is "correctly between" ref$AA_T_Mean and ref$AB_T_Mean
            else if (th <= tAB) {
                #interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
                BAF[i] <- 0.5* (th - tAA) / (tAB - tAA)
            }
            # 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be 
            # between ref$AB_T_Mean and ref$BB_T_Mean where it normally should not NOT BE.
            else {
                BAF[i] <- ifelse(th > med.tBB, 1, 0.5 + 0.5 * (th - tAB) / (med.tBB[i] - tAB))
            }
            
            if (th <= tAB) {
                eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
                LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
            }
            else {
                LRR[i] <- NaN
            }
        }
        # 7: positive for AA & BB & AB, Illumina style calculation
        else if (e.tAA & e.tAB & e.tBB){
            if (th < tAB) {
                BAF[i] <- ifelse(th < tAA, 0, 0.5 * (th - tAA) / (tAB - tAA))
                eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
            }
            else {
                BAF[i] <- ifelse(th >= tBB, 1, 0.5 + 0.5 * (th - tAB) / (tBB - tAB))
                eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
            }
            LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
        }
    }
    
    cbind(baf=.bound(BAF, 0, 1), lrr=LRR)
}

.theta <- function(data) 2/pi*atan(data$Y/data$X)

.bound <- function(x, low, high) {
    n <- .is.num(x)
    x[n & x < low] <- low
    x[n & x > high] <- high
    x
}

.my.mean <- function(x, minlen) {
    x <- x[.is.num(x)]
    l <- length(x)
    if (l < minlen) {
        NA
    }
    else if (l == 1) {
        x
    }
    else {
        mean(x)
    }
}

.my.sd <- function(x, minlen) {
    x <- x[.is.num(x)]
    l <- length(x)
    if (l < minlen) {
        NA
    }
    else if (l == 1) {
        0
    }
    else {
        sd(x)
    }
}

# Performs genomic wave correction using the same method as in
# PennCNV genomic_wave.pl --adjust mode. Bascially: fit a linear model of 
# LRR vs GC and then adjust the LRRs using the residuals.
# Reference: Diskin et al., Adjustment of genomic waves in signal intensities from 
#     whole-genome SNP genotyping platforms, Nucleic Acids Research 36:e126, 2008
.wave.correction <- function(metrics, gcmodel, distance=1000000) {
    # First train the model using only autosomal markers for which
    # LRR and GC fall within certain bounds (hardcoded in the perl script)
    gcmodel.auto <- gcmodel[gcmodel$Chr %in% 1:19,]
    w <- .is.num(metrics$lrr) & metrics$lrr > -1 & metrics$lrr < 1 & rownames(metrics) %in% gcmodel.auto$Name
    metrics.auto <- metrics[w,]
    gcmodel.auto <- gcmodel.auto[match(rownames(metrics.auto), gcmodel.auto$Name),]
    #w <- gcmodel.auto$GC > 15 & gcmodel.auto$GC < 80
    #metrics.auto <- metrics.auto[w,]
    #gcmodel.auto <- gcmodel.auto[w,]
    
    x <- NULL
    y <- NULL
    for (chrm in unique(gcmodel.auto$Chr)) {
        idx <- which(gcmodel.auto$Chr == chrm)
        pos <- gcmodel.auto[idx, "Position"]
        if (!is.null(distance)) {
            prev <- 0
            for (i in 1:length(idx)) {
                if (pos[i] - prev > distance) {
                    prev <- pos[i]
                }
                else {
                    idx[i] <- NA
                }
            }
            idx <- idx[!is.na(idx)]
        }
        x <- c(x, gcmodel.auto[idx, "GC"])
        y <- c(y, metrics.auto[idx, "lrr"])
    }
    w <- x > 15 & x < 80 & y > -1 & y < 1
    coeff <- lm(y~x, data.frame(x=x[w], y=y[w]))$coefficients
    
    # Now adjust LRR values for GC content
    m <- match(rownames(metrics), gcmodel$Name)
    w <- !is.na(m) & .is.num(metrics$lrr)
    metrics[w, "lrr"] <- metrics[w, "lrr"] - coeff[1] - coeff[2] * gcmodel[m[w], "GC"]
    metrics
}
