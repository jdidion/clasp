#' Plot BAF and LRR values for a single sample.
#'
#' \code{data} must be a data.frame with one row per marker and 
#' with at least two columns: chrm and pos (marker chromosome and position).
#' If, \code{plot.baf=TRUE}, there must be a \code{baf} column, and if
#' \code{plot.lrr=TRUE}, there must be a \code{lrr} column. If there is
#' also a \code{call} column, points in the BAF plot will be colored according
#' to genotype call.
#'
#' @param name sample name.
#' @param sample.id sampleID.
#' @param marker.set.id marker set ID.
#' @param metrics data.frame with metrics (BAF and LRR).
#' @param norm.params normalization params.
#' @param db database connection.
#' @param chromosomes the chromosomes to plot.
#' @param zoom if plotting a single chromosome, the region to zoom in on.
#' @param tick.interval distance between x-axis ticks (in Mb). Defaults to
#' a single tick per chromosome.
#' @param plot.baf whether to plot B allele frequencies.
#' @param plot.lrr whether to plot Log R Ratios.
#' @param lrr.ylim y-axis range within which to plot LRRs.
#' @param plot.trend whether to plot a trend line for LRR values.
#' @param density.band.window if non-null, overlays a plot of the region that contains
#' the specified fraction of values (\code{density.band.frac}) within windows the
#' specified size.
#' @param density.band.frac the fraction of values to include in the density band.
#' @param other.plots plotting function(s) to add additional panels to the BAF and LRR plots.
#' @param other.data additional data required by the additional plotting functions.
#' @param ... additional arguments to pass to panel plotting functions.
metrics.plot <- function(name, sample.id, marker.set.id, metrics, norm.params, db, chromosomes=1:19,
        zoom=NULL, tick.interval=NULL, plot.baf=TRUE, plot.lrr=TRUE, lrr.ylim=c(-2,2),
        plot.trend=TRUE, density.band.window=NULL, density.band.frac=0.75, 
        other.plots=NULL, other.data=NULL, ...) {

    calls <- dbGetPreparedQuery(db, .wrap("select marker_id, chromosome as chrm, position_bp as pos, call from 
        genotype join marker on genotype.marker_id = marker.id where genotype.sample_id = ?"),
        data.frame(sample.id))
    
    # HACK: see compute.metrics.Illumina for details.
    if ("switched" %in% names(norm.params) && any(norm.params$switched)) {
        w <- as.character(calls$marker_id) %in% names(norm.params$switched)[norm.params$switched]
        A <- w & calls$call == 1
        B <- w & calls$call == 2
        calls[A, "call"] <- 2
        calls[B, "call"] <- 1
    }

    # Only plot SNPs on specified chromosomes
    if (!is.null(chromosomes)) {
        w <- !is.na(match(calls$chrm, chromosomes))
        calls <- calls[w,]
        metrics <- metrics[w,]
    }
    chromosomes <- unique(calls$chrm)

    # Convert chromosome positions to absolute genomic positions
    chrms <- factor(calls$chrm, levels=chromosomes)
    chrm.sizes <- as.numeric(tapply(calls$pos, chrms, max))
    genome.size <- sum(chrm.sizes)
    names(chrm.sizes) <- chromosomes
    cs <- c(0, cumsum(as.numeric(chrm.sizes[-length(chrm.sizes)])))
    names(cs) <- chromosomes
    x <- cs[as.character(calls$chrm)] + calls$pos

    xlim <- c(0, genome.size)
    if (length(chromosomes) == 1 && !is.null(zoom)) {
        xlim <- zoom
    }

    nplots <- plot.baf + plot.lrr
    if (!is.null(other.plots)) {
        nplots <- nplots + length(other.plots)
    }
    if (nplots > 1) {
        par(mfrow=c(nplots, 1), mar=c(0.25,5,0.25,1), oma=c(3, 0, 3, 0))
    }

    if (plot.baf) {
        .plot.baf(x, calls, metrics$baf, chrms, xlim=xlim, ...)
    }

    if (plot.lrr) {
        .plot.scatter.with.trend.and.outliers(x, metrics$lrr, chrms, lrr.ylim, xlim=xlim,
            ylab="Log R Ratio", density.band.window, density.band.frac, plot.trend, ...)
    }
    
    if (!is.null(other.plots)) {
        for (i in 1:length(other.plots)) {
            f <- other.plots[[i]]
            o <- other.data[[i]]
            f(calls, metrics, o, chromosomes=cs, xlim=xlim, ...)
        }
    }

    if (is.null(tick.interval)) {
        tick.pos <- cs[1:length(chromosomes)]
        tick.name <- chromosomes
    }
    else {
        tick.pos <- NULL
        tick.name <- NULL
        for (i in 1:length(cs)) {
            ticks <- seq(1, chrm.sizes[i], tick.interval)
            tick.pos <- c(tick.pos, cs[i] + ticks)
            tick.name <- c(tick.name, chromosomes[i], ticks[2:length(ticks)] / 1000000)
        }
    }

    if (!is.null(name)) {
        title(name, outer=TRUE)
    }

    axis(1, tick.pos, tick.name)
}

.plot.baf <- function(x, calls, baf, chrms, A.col=c("lightblue","darkblue"), B.col=c("pink", "red"),
        H.col=c("mediumpurple1", "purple4"), N.col=c("lightgray", "darkgray"), mid.col='black', ...) {
    plot(0:1, type="n", ylim=c(0, 1), xaxt="n", xlab="", ylab="B Allele Frequency", ...)
    segments(0, 0.5, max(x), 0.5, col=mid.col, lty=2)

    odd <- as.integer(chrms) %% 2 == 1
    even <- as.integer(chrms) %% 2 == 0

    col <- rep(N.col[1], length(x))

    A <- calls$call == 1
    B <- calls$call == 2
    H <- calls$call == 3

    col[even] <- N.col[2]
    col[odd & A] <- A.col[1]
    col[even & A] <- A.col[2]
    col[odd & H] <- H.col[1]
    col[even & H] <- H.col[2]
    col[odd & B] <- B.col[1]
    col[even & B] <- B.col[2]

    segments(x[1], 0.5, x[length(x)], 0.5, lty=2)
    points(x, baf, col=col, pch=20)
}

.plot.scatter.with.trend.and.outliers <- function(x, y, chrms, ylim, density.band.window,
        density.band.frac, plot.trend=!is.null(density.band.window),
        point.col=c("lightgray", "darkgray"), zero.col="black",
        trend.col="red", density.fill.col=rgb(1,0,0,0.33), ...) {
    odd <- as.integer(chrms) %% 2 == 0
    even <- as.integer(chrms) %% 2 == 1

    outliers.low <- y < ylim[1]
    y[outliers.low] <- ylim[1]
    outliers.high <- y > ylim[2]
    y[outliers.high] <- ylim[2]

    col <- rep(point.col[1], length(x))
    col[even] <- point.col[2]
    col[outliers.high | outliers.low] <- trend.col

    plot(0:1, type="n", ylim=ylim, xaxt="n", ...)
    points(x, y, col=col, pch=20)

    segments(0, 0, max(x), 0, col=zero.col, lty=2)

    if (plot.trend) {
        # Plot "trend lines" for each chromosome. Currently using the supersmoother function,
        # but could also use a kernel smoother, regression, splines or loess. See:
        # http://statistics.berkeley.edu/classes/s133/Smooth-a.html
        # http://statistics.berkeley.edu/classes/s133/Smooth1a.html
        fits <- by(data.frame(x=x, y=y), chrms, function(m) {
          w <- .is.num(m$x) & .is.num(m$y)
          supsmu(m$x[w], m$y[w])
        })
        for (n in names(fits)) {
            lines(fits[[n]]$x, fits[[n]]$y, col=trend.col, lwd=2)
        }
    }
    if (!is.null(density.band.window)) {
        # Use a sliding window to determine the upper and lower bounds of y for density.band.frac SNPs
        bounds <- by(data.frame(x=x, y=y), chrms, function(m) {
            # Calculate the upper and lower quantiles for all windows
            rem <- (1 - density.band.frac) / 2
            if (nrow(m) <= density.band.window) {
                q <- quantile(m$y, c(rem, 1 - rem), na.rm=T)
                m$low <- q[1]
                m$high <- q[2]
            }
            else {
                last.win <- nrow(m) - density.band.window + 1
                m$low <- NA
                m$high <- NA
                m[1:last.win, c("low", "high")] <- t(sapply(seq(1, last.win), function(start) {
                    end <- start + density.band.window - 1
                    y <- m$y[start:end]
                    y <- y[.is.num(y)]
                    quantile(y, c(rem, 1 - rem))
                }))
                # Copy the value of the last window to cover the remaining x values
                m[(last.win + 1):nrow(m), c("low", "high")] <- m[last.win, c("low", "high")]
            }
            m
        }, simplify=FALSE)

        # Draw filled band lines
        for (n in names(bounds)) {
            m <- bounds[[n]]
            lines(m$x, m$low, col=trend.col, lwd=0.5)
            lines(m$x, m$high, col=trend.col, lwd=0.5)
            polygon(c(m$x, rev(m$x)), c(m$low, rev(m$high)), border=NA, col=density.fill.col)
        }
    }
}

#' Apply an alpha adjustment to a color
#'
#' @param r color name (convertible by col2rgb), or vector or matrics of
#' RGB values.
#' @param alpha alpha value to apply.
#' 
#' @return vector of hex colors.
rgb.alpha <- function(r, alpha=1.0) {
    if (is.character(r)) {
        r <- col2rgb(r)
    }
    if (is.matrix(r)) {
        if (length(alpha) == 1) {
            apply(r, 2, rgb.alpha, alpha)
        }
        else {
            sapply(1:ncol(r), function(i) rgb.alpha(r[,i], alpha[i]))
        }
    }
    else {
        rgb(r[1], r[2], r[3], alpha * 255, maxColorValue=255)
    }
}