# Recode genotypes as bit vectors (integers).
#
# Character genotypes are converted to bit vectors. Currently, CLASP only supports
# bi-allelic markers, so the bit vectors are of length two, i.e. the possible values are:
# 00=0: missing information
# 01=1: reference allele
# 10=2: variant allele
# 11=3: heterozygous
# 
# @param geno A vector of character genotypes.
# @param alleles A two-element vector: the reference and variant alleles.
# @param N.call A vector of the possible character value(s) of missing information.
#
# @return
# A vector of recoded genotypes.
.recode <- function(geno, alleles, N.calls="N") {
    recoded <- rep(0, length(geno))
    recoded[geno == alleles[1]] <- 1
    recoded[geno == alleles[2]] <- 2
    recoded[recoded == 0 & !(geno %in% N.calls)] <- 3
    recoded
}

# Determine which markers have consistent genotypes across a set of replicate samples.
#
# @param geno A MxN matrix of (recoded) genotypes, where M = markers and N = samples.
#
# @return
# A logical vector of length M.
.consistent <- function(geno, ignore.missing=TRUE, allow.hets=TRUE) {
    apply(geno, 1, function(x) {
        if (ignore.missing) {
            x = x[x>0]
        }
        if (allow.hets) {
            x <- x[x<3]
        }
        length(unique(x)) <= 1
    })
}

# Compute allele and genotype frequencies. These are required as parameters to several
# other functions.
#
# @param g A vector of (recoded) genotypes.
#
# @return
# A list with the following elements:
# AA, BB, AB, N: number of homozygous reference, homozygous variant, heterozygous and no-call
# genotypes, respectively.
# A, B: number of reference and variant alleles, respectively.
# A.counts, B.counts: number of A and B alleles, respectively, in each genotype, i.e. integer
# vectors of the same length as \code{g} and possible values 0,1,2.
# p, q, N.freq: A, B and N allele frequencies.
# n.geno, n.alleles: Total number of genotypes (i.e. \code{length(g)}) and total number of (non-N)
# alleles, respectively.
.allele.freq <- function(g) {
    AA <- g == 1
    BB <- g == 2
    AB <- g == 3
    A <- (2 * AA) + AB
    B <- (2 * BB) + AB
    a <- sum(A)
    b <- sum(B)
    tot <- a + b
    N <- sum(g == 0)
    n <- length(g)
    list(AA=sum(AA), BB=sum(BB), AB=sum(AB), A=a, B=b, N=N, N.freq=N/n,
         A.counts=A, B.counts=B, n.geno=n, n.alleles=tot, p=a/tot, q=b/tot)
}

# Compute expected allele frequencies and do a Pearson's Chi Squared test with HWE as the
# null hypothesis.
#
# @param af A list of lists, where each element is the result of a call to \code{allele.freq}.
#
# @return A list of two elements:
# 1. x2: chi-squared value (one df)
# 2. pvalue: p-value of the Chi Squared test.
.hwe.test <- function(af) {
    EAA <- (af$p ^ 2) * af$n.alleles
    EBB <- (af$q ^ 2) * af$n.alleles
    EAB <- 2 * af$p * af$q * af$n.alleles
    x2 <- ((af$AA - EAA) / EAA) + ((af$BB - EBB) / EBB) + ((af$AB - EAB) / EAB)
    list(x2=x2, pvalue=dchisq(x2,1))
}

.pairwise.distances <- function(geno) {
	nsamples <- ncol(geno)
    pairs <- combn(1:nsamples, 2)
    dists <- sapply(1:ncol(pairs), function(i) sum(geno[,pairs[1,i]] != geno[,pairs[2,i]]))
	data.frame(t(pairs), dists)
}

.distance.matrix <- function(dists, sample.names) {
	n <- length(sample.names)
	distmat <- matrix(0, n, n, dimnames=list(sample.names, sample.names))
    for (i in 1:nrow(dists)) {
        distmat[dists[i,1],dists[i,2]] <- dists[i,3]
        distmat[dists[i,2],dists[i,1]] <- dists[i,3]
    }
	distmat
}

.orthoganal.distance.matrix <- function(g1, g2, dnames) {
	odists <- do.call(rbind, lapply(1:ncol(g1), function(i)
        sapply(1:ncol(g2), function(j) sum(g1[,i] != g2[,j]))))
  dimnames(odists) <- dnames
	odists
}

# Compute the pairwise LD (\code{r^2} value) between all pairs of markers.
#
# @param af Allele frequency information for the set of markers. A list of lists, where each
# element is the result of a call to \code{allele.freq}.
.LD <- function(af) {
    pairs <- t(data.frame(combn(1:length(af), 2)))
    rownames(pairs) <- 1:nrow(pairs)
    ld <- do.call(rbind, apply(pairs, 1, function(pair) {
        g1 <- af[[pair[1]]]
        g2 <- af[[pair[2]]]

        p1 <- max(g1$p, g1$q)
        p2 <- max(g2$p, g2$q)
        q1 <- 1-p1
        q2 <- 1-p2

        maj1 <- which(c(g1$p, g1$q) == p1)[1]
        maj2 <- which(c(g2$p, g2$q) == p2)[1]

        Dmin <- max(-p1*p2, -q1*q2)
        pmin <- p1*p2 + Dmin;
        Dmax <- min(p1*q2, p2*q1);
        pmax <- p1*p2 + Dmax;

        counts <- table(g1[[c("A.counts","B.counts")[maj1]]], g2[[c("A.counts","B.counts")[maj2]]])

        n3x3 <- matrix(0, nrow=3, ncol=3)
        colnames(n3x3) <- rownames(n3x3) <- 0:2
        # ensure the matrix has highest frequency values in upper left
        for (i in rownames(counts))
            for (j in colnames(counts))
                n3x3[3-as.numeric(i), 3-as.numeric(j)] <- counts[i,j]

        loglik <- function(pAB, ...) {
            (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
            (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(p1-pAB) +
            (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(p2-pAB) +
            (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-p1-p2+pAB) +
            n3x3[2,2]*log(pAB*(1-p1-p2+pAB) + (p1-pAB)*(p2-pAB))
        }

        eps <- .Machine$double.eps
        pAB <- optimize(loglik, lower=pmin+eps, upper=pmax-eps, maximum=TRUE)$maximum

        estD <- pAB - p1*p2
        if (estD > 0)
            estDp <- estD / Dmax
        else
            estDp <- estD / Dmin

        n <-  sum(n3x3)
        m <- p1 * q1 * p2 * q2
        corr <- estD / sqrt(m)
        dchi <- (2 * n * estD^2) / m
        dpval <- 1 - pchisq(dchi,1)

        data.frame(r2=corr^2, n=n, x2=dchi, pvalue=dpval)
    }))

    list(pairs=pairs, ld=ld)
}

# Generate genotypes for a synthetic cross between two samples.
#
# \code{genotypes} is a matrix of at least two columns, where row names are marker IDs and
# column names are sample IDs. \code{sample.pairs} is either a list of vectors or a matrix that
# specifies the mother and father (in order) for each pair. If there is a third element, it is
# used as the name for the synthetic F1. \code{genotypes} only has two columns and sample.pairs
# is NULL, the first column is assumed to be the mother and the second the father.
#
# @param genotypes A MxN matrix specifying all parental genotypes. See details.
# @param sample.pairs A list or matrix specifying the mother and father for each cross.
.impute.f1.genotypes <- function(genotypes, sample.pairs=NULL) {
    if (is.null(sample.pairs)) {
        sample.pairs <- t(combn(colnames(genotypes), 2))
    }
    else if (!(is.matrix(sample.pairs) | is.data.frame(sample.pairs))) {
        sample.pairs <- do.call(rbind, sample.pairs)
    }

    # TODO: this is pretty fast, but we could parallelize it
    synth <- apply(sample.pairs[,1:2], 1, function(x)
        .impute.f1(genotypes[,as.character(x[1])], genotypes[,as.character(x[2])]))
    rownames(synth) <- rownames(genotypes)
    if (ncol(sample.pairs) > 2) {
        colnames(synth) <- sample.pairs[,3]
    }
    list(sample.pairs=sample.pairs, genotypes=synth)
}

.impute.f1 <- function(a, b) {
    newg <- rep(0, length(a))

    aA <- a==1
    aB <- a==2
    aH <- a==3
    bA <- b==1
    bB <- b==2
    bH <- b==3
    AH <- (aA & bH)|(bA & aH)
    BH <- (aB & bH)|(bB & aH)
    HH <- aH & bH

    newg[aA & bA] <- 1
    newg[aB & bB] <- 2
    newg[(aA & bB)|(aB & bA)] <- 3
    newg[AH] <- sample(c(1,3), sum(AH), TRUE)
    newg[BH] <- sample(c(2,3), sum(BH), TRUE)
    newg[HH] <- sample(c(1,2,3,3), sum(HH), TRUE)

    newg
}

.f1.consistent <- function(geno, ignore.missing=TRUE) {
    apply(geno, 1, function(x) {
        if (any(x == 0)) {
            ignore.missing || all(x == 0) || (sum(x[2:3] == 0) == 1)
        }
        else if (all(x == 1) || all(x == 2)) {
            TRUE
        }
        else if (x[1] == 3 && (x[2] == 3 && x[3] == 3 || x[2] != x[3])) {
            TRUE
        }
        else if ((x[1] == x[2] && x[3] == 3) || (x[1] == x[3] && x[2] == 3)) {
            TRUE
        }
        else {
            FALSE
        }
    })
}
