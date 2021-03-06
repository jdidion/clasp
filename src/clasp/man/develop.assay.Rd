% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/assay.R
\name{develop.assay}
\alias{develop.assay}
\title{Develop a new assay.}
\usage{
develop.assay(name, marker.set, db, sample.set = NULL, insert = TRUE,
  description = NULL, ob.name = NULL, ob.description = NULL,
  exclude.cell.lines = TRUE, min.call.rate = 0.8, min.maf = 0,
  max.hwe.pvalue = NULL, consistency.check = TRUE, max.ld.r2 = NULL,
  max.ld.pvalue = NULL, keep.informative.linked.markers = TRUE,
  max.sdp.diff = 0.05, diag.coeffs = list(N = 1, AA = 1, BB = 1, AB = 1),
  reduce = FALSE, min.diff = NULL, min.markers = NULL,
  diag.threshold.iters = seq(0.25, 0.7, 0.05))
}
\arguments{
\item{name}{Assay name}

\item{marker.set}{A vector of marker set IDs (see details) or a marker data.frame.}

\item{db}{The database connection.}

\item{sample.set}{A vector of sample set IDs, a sample data.frame or NULL if all samples should be used.}

\item{insert}{Whether or not to insert the new assay into the database.}

\item{description}{If insert=TRUE, the text description of the new marker set.}

\item{ob.name}{If insert=TRUE, the name of the new outbred-specific marker set.}

\item{ob.description}{If insert=TRUE, the text description of the new outbred-specific marker set.}

\item{exclude.cell.lines}{If TRUE, cell line samples are not used to develop the assay.}

\item{min.call.rate}{Minimum fraction of samples that have to have a non-missing (i.e. non-zero)
genotype (inclusive).}

\item{min.maf}{Minimum fraction of samples with non-missing genotype information that must have
a minor allele. The value is not inclusive, so a value of 0 means the MAF must be > 0, whereas a
value of NULL means the check is not performed.}

\item{max.hwe.pvalue}{The maximum p-value (inclusive) for a test of deviation from the expected
Hardy-Weinberg allele frequency.}

\item{consistency.check}{Whether a consistency check should be performed for replicates.}

\item{max.ld.r2}{Maximum \code{R^2} value (inclusive) for pairwise linkage between consecutive markers.}

\item{max.ld.pvalue}{Maximum p-value (inclusive) for a test of pairwise linkage between
consecutive markers.}

\item{keep.informative.linked.markers}{If a linkage test is performed, this parameter specifies
that markers that exceed the \code{R^2} threshold should not be filtered if their SDPs differ
(see Details).}

\item{max.sdp.diff}{Maximum fraction by which two SDPs can differ and still be considered equal
for the purposes of \code{keep.informative.linked.markers}.}

\item{diag.coeffs}{genotype-speific coefficients for computing diagnostic values.}

\item{reduce}{Logical indicating whether the maximal marker set should be pruned using the
the following parameters.}

\item{min.diff}{When reduce=TRUE, the minimum allowable distance between any two markers.}

\item{min.markers}{When reduce=TRUE, the minimum number of markers in the final marker set. The
actual size of the final marker set may be smaller than this if the other filtering parameters
remove too many markers.}

\item{diag.threshold.iters}{A vector of thresholds to iterate through when pruining markers.}
}
\value{
A list with the following elements:
set.id: The ID of the new marker set created with the markers for this assay, or NULL if insert=FALSE.
markers: Full information for the set of markers included in the assay.
samples: Full information for the set of all samples considered for this assay.
sample.ixs: Indices of samples that were actually used to develop the assay (i.e. after removing
replicates).
replicates: A list  of replicate groups.
geno: MxN matrix of genotypes, where M is the assay markers and N is the assay samples.
diagnostic.values: MxN matrix of the diagnostic value for each SNP and for each sample.
distmat: NxN matrix of the pairwise number of unequal markers.
}
\description{
Develop a new assay by selecting the most informative subset of \code{markers} for the
specified samples.
}
\details{
An assay may be developed from a single marker set or from multiple intersecting marker sets. In
the latter case (i.e. if \code{markers} is a vector of two or more IDs), the intersection of all
the marker sets will be found using \code{intersect.marker.sets}.

For the purposes of the consistency check and removing replicate samples, two samples are
considered replicates if they meet any of the following criteria:
1. Identical names
2. Identical non-NULL cell line names AND identical sources
3. Null cell line names AND identical non-NULL backgrounds

An SDP (sample distribution pattern) is the partitioning of a set of samples based on genotypes
for a single marker. It can be thought of as the phylogenetic tree with N braches, where N is the
number of possible genotypes. Consecutive markers with different SDPs are considered informative
even if their degree of linkage exceeds the specified \code{R^2} threshold.

It is expected that all samples in the database have passed QC checks; this function does not
perform any sample-level QC.
}

