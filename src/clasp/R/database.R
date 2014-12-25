#' Open the SQLite database specified by \code{dbfile}.
#'
#' @param dbfile the database file
#'
#' @return a datbase connection
connect.database <- function(dbfile) {
    dbConnect(SQLite(), dbname=dbfile)
}

#' Initialize a new database with the default schema.
#'
#' Using an alternate schema definition will probably break the program
#' unless you make appropriate changes to R the code as well.
#'
#' @param dbfile the database file to create
#' @param schema.file the file from which to load the DDL.
#'
#' @return The open database connection.
init.database <- function(dbfile, schema.file=NULL) {
    if (is.null(schema.file)) {
        schema.file <- system.file("schema.SQL", package="clasp")
    }
    db <- connect.database(dbfile)
    ddl <- strsplit(paste(readLines(schema.file), collapse="\n"), ";", fixed=TRUE)[[1]]
    for (stmt in ddl) {
        dbSendQuery(db, stmt)
    }
    invisible(db)
}

#' Define a genome build in the database.
#'
#' For organisms that have multiple maintainers of the same build, the most official name
#' should be used. E.g. for mouse, 'GRCm38' should be used rather than 'mm10'.
#'
#' @param taxon Taxon of the genome (i.e. Mus musculus)
#' @param version Genome build version.
#' @param db database connection
#'
#' @return the ID of the genome.
import.genome <- function(taxon, version, db) {
    tryCatch({
        dbBegin(db)
        dbSendPreparedQuery(db, "INSERT OR IGNORE INTO Genome (taxon, version) VALUES (?,?)",
            data.frame(taxon, version))
        genome.id <- dbGetPreparedQuery(db, "SELECT id FROM Genome WHERE taxon = ? AND version = ?",
            data.frame(taxon, version))[1,1]
        dbCommit(db)
        return(genome.id)
    }, error=function(e) {
        dbRollback(db)
        stop(e)
    })
}

#' Get information for all genome builds in the database.
#'
#' @param db The database connection.
#'
#' @return table of genome information.
get.all.genomes <- function(db) {
    dbGetQuery("SELECT * FROM Genome")
}

#' Get genome ID by name.
#'
#' @param name genome name
#' @param db database connection
#'
#' @return the ID
get.genome.id <- function(name, db) {
    dbGetQuery(db, "SELECT id FROM Genome WHERE name = ?", name)
}

#' Define a genotyping array in the database (if it
#' doesn't already exist) and import the default marker 
#' set for that array for a given genome.
#'
#' @param name array name
#' @param platform array platform (e.g. Illumina)
#' @param genome.id ID of the genome on which the specified
#' marker annotations are based.
#' @param markers marker data, as described in import.marker.set
#' @param db the database
#' @param ... additional parameters to import.marker.set
#'
#' @return a vector of two IDs: (array.id, marker.set.id)
import.array <- function(name, platform, genome.id, markers, db, ...) {
    tryCatch({
        dbBegin(db)
        dbSendPreparedQuery(db, "INSERT OR IGNORE INTO Array (name, platform) VALUES (?,?)",
            data.frame(name, platform))
        array.id <- dbGetPreparedQuery(db, "SELECT id FROM Array WHERE name = ? AND platform = ?",
            data.frame(name, platform))[1,1]
        marker.set <- import.marker.set(name, markers, db, genome.id=genome.id, start.transaction=FALSE, ...)
        dbSendPreparedQuery(db, "INSERT INTO ArrayXGenome VALUES(?,?,?)", 
            data.frame(array.id, genome.id, marker.set$set.id))
        dbCommit(db)
        return(list(array.id=array.id, marker.set=marker.set))
    }, error=function(e) {
        dbRollback(db)
        stop(e)
    })
}

#' Get information for all genotyping arrays in the database.
#'
#' @param db the database connection.
#'
#' @return table of array information.
get.all.arrays <- function(db) {
    dbGetQuery("SELECT * FROM Array")
}

#' Get array ID by name.
#'
#' @param name the array name.
#' @param db database connection
#'
#' @return the ID.
get.array.id <- function(name, db) {
    dbGetQuery(db, "SELECT id FROM Array WHERE name = ?", name)
}

#' Import marker data.
#'
#' Currently, CLASP only supports bi-allelic markers and does not support sex chromosome markers.
#' The later is because cell lines are often of indeterminate sex or exhibit sex chromosome mosaicism.
#'
#' If \code{markers} is a vector, the elements are marker IDs and the names can optionally be marker
#' set-specific names. If it is a data.frame, it must have five columns:
#' 1. name: Name of the marker.
#' 2. chromosome: Name of the chromosome. Must be an integer.
#' 3. bp: Physical position of the marker.
#' 4. cm: Genetic position of the marker (in centimorgans).
#' 5. alleles: String listing the two possible genotypes for this marker. The alleles should be listed
#' with the reference allele first, followed by the variant allele. If this marker is to be
#' associated with intensity data, the intensity values should be reordered to match the allele order.
#'
#' @param name Name for the marker set.
#' @param genome.id ID of the genome build from which the markers were sourced.
#' @param markers vector of marker IDs or path to marker data file or data.frame containing marker 
#' information. See Details.
#' @param db the datbase connection.
#' @param description Text description of the marker set.
#' @param start.transaction set to false if calling from another method that starts a transaction.
#'
#' @return
#' If markers is a path, returns a list with two elements, set.id and markers. Otherwise returns the 
#' integer ID of the new marker set.
import.marker.set <- function(name, markers, db, description=NA, genome.id=NULL, start.transaction=TRUE) {
    tryCatch({
        if (start.transaction) dbBegin(db)
        dbSendPreparedQuery(db, .wrap("INSERT INTO MarkerSet (name, description) VALUES (?,?)"),
            data.frame(name, description))
        set.id <- dbGetQuery(db, "SELECT last_insert_rowid() FROM MarkerSet")[1,1]
		result <- set.id
        if (!is.data.frame(markers)) {
	        if (length(markers) == 1) {
	        	markers <- read.table(markers, sep="\t", header=T, stringsAsFactors=F)
	        }
			else {
				if (is.null(names(markers))) {
					names(markers) <- paste(name, seq(to=length(markers)), sep=".")
				}
	            dbSendPreparedQuery(db, "INSERT INTO MarkerXMarkerSet VALUES (?,?,?)", data.frame(
	                markers, rep(set.id, length(markers)), names(markers), stringsAsFactors=FALSE))
	        }
        }
		if (is.data.frame(markers)) {
            if (is.character(genome.id)) {
                genome.id <- get.genome.id(genome.id)
            }
			marker.ids <- rep(NA, nrow(markers))
            for (i in 1:nrow(markers)) {
                dbSendPreparedQuery(db, .wrap("INSERT OR IGNORE INTO Marker (genome_id,
                    chromosome, position_bp, position_cm, alleles) VALUES (?,?,?,?,?)"),
                    as.data.frame(cbind(genome.id, markers[i,2:5])))
                marker <- dbGetPreparedQuery(db, .wrap("SELECT id, alleles FROM Marker WHERE
                    genome_id = ? AND chromosome = ? AND position_bp = ?"),
                    as.data.frame(cbind(genome.id, markers[i,2:3])))
                if (nrow(marker) > 1 || marker[1,2] != markers[i,5]) {
                    stop(paste("Duplicate marker with different alleles", marker[1,1]))
                }
				marker.ids[i] <- marker[1,1]
                dbSendPreparedQuery(db, "INSERT INTO MarkerXMarkerSet VALUES (?,?,?)",
                    data.frame(marker.ids[i], set.id, markers[i,1]))
            }
			result <- list(set.id=set.id, markers=cbind(id=marker.ids, markers))
        }
        if (start.transaction) dbCommit(db)
        return(result)
    }, error=function(e) {
        if (start.transaction) dbRollback(db)
        stop(e)
    })
}

#' Find the intersection of two or more marker sets.
#'
#' @param marker.set.ids A vector of marker.set.ids.
#' @param db The database connection.
#' @param insert Logical, whether the intersected marker set should be inserted into the database.
#' @param name Name for the new marker set.
#' @param description Text description of the new marker set.
#'
#' @return
#' A data frame with 4+N columns, where N is the number of marker sets. The first four columns are:
#' chromosome, position_bp, position_cm, alleles. The remaining columns are the marker set-specific
#' names for each marker. If insert=TRUE, the result is actually a two-element list: set.id=new marker 
#' set id, markers=marker data frame.
intersect.marker.sets <- function(marker.set.ids, db, insert=FALSE, name=NULL, description=NULL) {
    sql <- paste(.wrap("SELECT Marker.*, GROUP_CONCAT(MarkerXMarkerSet.name, ','),
        COUNT(DISTINCT MarkerXMarkerSet.set_id) AS cnt FROM Marker JOIN MarkerXMarkerSet
        ON Marker.id = MarkerXMarkerSet.marker_id WHERE MarkerXMarkerSet.set_id IN ("),
        paste(marker.set.ids, collapse=','), ") GROUP BY Marker.id HAVING cnt =",
        length(marker.set.ids), "ORDER BY chromosome, position_bp")
    markers <- dbGetQuery(db, sql)
    markers <- data.frame(markers[,c(1,3,4,5,6)], do.call(rbind,
        strsplit(markers[,7:(ncol(markers)-1)], ",", fixed=TRUE)), stringsAsFactors=FALSE)
    names <- dbGetQuery(db, paste("SELECT id,name FROM MarkerSet WHERE id IN (",
        paste(marker.set.ids, collapse=','), ")"))
    colnames(markers)[6:ncol(markers)] <- names[match(marker.set.ids, names[,1]),2]
    if (insert) {
        tryCatch({
            dbBegin(db)
            dbSendPreparedQuery(db, "INSERT INTO MarkerSet (name, description) VALUES (?,?)",
                data.frame(name, description))
            set.id <- dbGetQuery(db, "SELECT last_insert_rowid() FROM MarkerSet")[1,1]
            dbSendPreparedQuery(db, "INSERT INTO MarkerXMarkerSet (marker_id, set_id) VALUES (?,?)",
                data.frame(markers[,1], rep(set.id, nrow(markers))))
            dbCommit(db)
        }, error=function(e) {
            dbRollback(db)
            stop(e)
        })
        markers <- list(set.id=set.id, markers=markers)
    }
    invisible(markers)
}

#' Get a list of all marker sets.
#'
#' @param db the datbase connection.
#'
#' @return a table of marker sets.
get.all.marker.sets <- function(db) {
    dbGetQuery(db, "SELECT * FROM MarkerSet")
}

#' Get a specific marker set.
#'
#' @param set.id The marker set ID.
#' @param db the datbase connection.
#'
#' @return table with marker information.
get.marker.set <- function(set.id, db) {
    dbGetPreparedQuery(db, "SELECT * FROM MarkerSet WHERE id = ?", data.frame(set.id))
}

#' Load the markers from a marker set or list of marker IDs.
#'
#' @param ids The marker set ID, or a list of marker IDs.
#' @param db the datbase connection.
#'
# @return a data frame with 5 columns: chromosome, position_bp and position_cm, alleles, 
#' name (the marker set-specific name for the marker, if any).
get.markers <- function(ids, db) {
	if (length(ids) == 1) {
	    markers <- dbGetPreparedQuery(db, .wrap("SELECT id,chromosome,position_bp,position_cm,alleles,
	        MarkerXMarkerSet.name FROM Marker JOIN MarkerXMarkerSet ON Marker.id =
	        MarkerXMarkerSet.marker_id WHERE MarkerXMarkerSet.set_id = ?"), data.frame(ids))
	}
	else {
		markers <- dbGetPreparedQuery(db, .wrap("SELECT id,chromosome,position_bp,position_cm,alleles
			FROM Marker WHERE id = ?"), data.frame(ids))
	}
    invisible(markers)
}

#' Create a sample set from a list of samples.
#'
#' A sample set is created. If samples is a vector of sample IDs, those samples are added to
#' the set. If it is a data frame, \code{import.samples} is called.
#'
#' @param name A unique name for the sample set.
#' @param db the datbase connection.
#' @param description A text description of the sample set.
#' @param samples a vector of sample IDs, a data.frame as specified in \code{import.samples},
#' or a path to a file containing sample data.
#' @param array.id array ID to use for import.samples. Ignored unless samples is a data frame
#' or path.
#' @param pedigree Pedigree information as specified in \code{import.samples}.
#' @param fetch if \code{samples} is a vector of sample IDs, whether to fetch the sample information
#' from the database rather than the default, which is to just return the new sample set ID.
#'
#' @return sample set ID, or, if samples was a data frame, a list with the following elements:
#' 1. set.id: The ID of the new sample set.
#' 2. ids: A vector of ids for the samples that were imported.
#' 3. ped: Pedigree data.frame, with sample names replaced with IDs.
#'
#' @seealso \code{\link{import.samples}}
import.sample.set <- function(name, db, description=NULL, samples=NULL, array.id=NULL, pedigree=NULL, fetch=FALSE) {
    tryCatch({
        dbBegin(db)
        dbSendPreparedQuery(db, "INSERT INTO SampleSet (name, description) VALUES (?,?)",
            data.frame(name, description))
        set.id <- dbGetQuery(db, "SELECT last_insert_rowid() FROM SampleSet")[1,1]
        result <- NULL
        if (!is.null(samples)) {
			if (is.data.frame(samples) || length(samples)==1) {
				result <- list(set.id=set.id, samples=import.samples(samples, array.id, db, pedigree, start.transaction=FALSE))
                sample.ids <- result$samples$id
            }
            else {
                sample.ids <- samples
            }
            dbSendPreparedQuery(db, "INSERT INTO SampleXSampleSet VALUES (?,?)",
                data.frame(sample.ids, rep(set.id, length(sample.ids))))
        }
        if (is.null(result)) {
            if (fetch) {
                result <- list(set.id=set.id, samples=get.samples(set.id, db))
            }
            else {
                result <- set.id
            }
        }
        dbCommit(db)
        invisible(result)
    }, error=function(e) {
        dbRollback(db)
        stop(e)
    })
}

#' Import sample data and add the samples to an existing sample set.
#'
#' @param sample.set.id The sample set to add to.
#' @param samples a data.frame as specified in \code{import.samples}.
#' @param array.id ID of the array.
#' @param db The database connection.
#' @param pedigree Pedigree information as specified in \code{import.samples}.
#'
#' @return a list with two elements:
#' 1. ids: a vector of database IDs, one for each sample.
#' 2. ped: the pedigree matrix with sample names replaced by ids (or NULL if no pedigree was specified).
add.samples <- function(sample.set.id, samples, array.id, db, pedigree=NULL) {
    result <- import.samples(samples, array.id, db, pedigree)
    sample.ids <- result$samples$id
    tryCatch({
        dbBegin(db)
        dbSendPreparedQuery(db, "INSERT INTO SampleXSampleSet VALUES (?,?)",
            data.frame(sample.ids, rep(sample.set.id, length(sample.ids))))
        dbCommit(db)
    }, error=function(e) {
        dbRollback(db)
        stop(e)
    })
    invisible(result)
}

#' Import sample data.
#'
#' The samples and pedigree parameters can be specified as a single or two separate
#' data.frames or files. If they are specfied as files, than blank fields in logical
#' columns will be given the default value of FALSE.
#'
#' \code{samples} must be a data.frame with eight columns:
#' 1. name: Sample name.
#' 2. cell_line: Name of the cell line, or blank if this is a non-cell line sample.
#' 3. background: Genetic background information, e.g. name of mouse inbred strain.
#' 4. type: sample type (inbred, outbred, F1, wild); see note.
#' 5. taxon: The sample taxon, as precise as possible.
#' 6. source: Person/lab/institution that provided the sample.
#' 7. description: Text description of the sample.
#' 8. impute: boolean, whether the sample should be used for imputation
#' of synthetic F1s.
#'
#' \code{pedigree} must be a data.frame with three columns: name, mother, father. These must all
#' uniquely refer to a sample within the \code{samples} data frame.
#'
#' If \code{samples} is a path, then pedigree can be specified as an additional three columns.
#'
#' Samples must be coded as being of one of four types:
#' 1. inbred: fully inbred mouse strain with little to no variablility between the animals of
#' that strain.
#' 2. outbred: outbred stock; individuals share a common background but are randomly mated such
#' that there is a large degree of variability between individuals.
#' 3. F1: a cross between two distinct parental strains/individuals. Pedigree information should be
#' provided if the parental strains are also among the samples.
#' 4. wild: wild-caught individual; generally does not exhibit high identity with any other strain
#' in the database.
#'
#' @param samples path to file, data.frame containing sample information, or vector of sample IDs.
#' @param array.id ID of the array.
#' @param db the database connection.
#' @param pedigree pedigree information for samples with known parentage, or NULL.
#' @param fetch if \code{samples} us a list of IDs, whether to fetch sample information.
#' @param start.transaction whether to start a database transaction. Should only be set to FALSE
#' when called from another function that itself starts a transaction.
#'
#' @return a list with two elements:
#' 1. ids: a vector of database IDs, one for each sample.
#' 2. ped: the pedigree matrix with sample names replaced by ids (or NULL if no pedigree was specified).
import.samples <- function(samples, array.id, db, pedigree=NULL, fetch=FALSE, start.transaction=TRUE) {
    if (is.null(array.id)) {
        stop("Must specify non-NULL array.id")
    }
    else if (is.character(array.id)) {
        array.id <- get.array.id(array.id)
    }
    tryCatch({
        if (start.transaction) dbBegin(db)
        now <- Sys.time()
        if (!is.data.frame(samples)) {
        	samples <- read.table(samples, sep="\t", header=T, stringsAsFactors=F, comment.char="",
				colClasses=c("character","character","character","character","character",
							 "character","character","logical","character","character"))
		}
		if (ncol(samples) > 8) {
			pedigree <- samples[,c(1,9,10)]
			samples <- samples[,1:8]
		}
        samples <- data.frame(.convert.samples(samples), array_id=array.id, time_added=rep(now, nrow(samples)))
        dbSendPreparedQuery(db, .wrap(
            "INSERT INTO Sample (name, cell_line, background, type, taxon, source, description, 
			impute, array_id, time_added) VALUES (?,?,?,?,?,?,?,?,?,?)"), samples)
        ids <- dbGetPreparedQuery(db, "SELECT id FROM Sample WHERE time_added = ? ORDER BY id ASC", 
            data.frame(now))[,1]
		samples <- data.frame(id=ids, samples)
        if (!is.null(pedigree)) {
            samples$mother <- NA
            samples$father <- NA
            # pedigree is provided as a mapping between sample names
            if (ncol(pedigree) == 3) {
                w <- nchar(pedigree[,2]) > 0 & nchar(pedigree[,3]) > 0
                ped <- NULL
                for (i in which(w)) {
                    s <- which(pedigree[i,1] == samples[,2])
                    m <- which(pedigree[i,2] == samples[,2])
                    f <- which(pedigree[i,3] == samples[,2])
                    if (length(s) != 1 || length(m) != 1 || length(f) != 1) {
                        stop(paste("Missing or non-unique sample", pedigree[i,]))
                    }
                    ped <- rbind(ped, c(ids[s], ids[m], ids[f]))
                }
                samples[w, c("mother", "father")] <- ped
            }
            # pedigree is provided as a 2xN matrix of mother/father IDs
            else if (ncol(pedigree) == 2 && nrow(pedigree) == nrow(samples)) {
                ped <- cbind(ids, pedigree)
                samples[, c("mother", "father")] <- pedigree
            }
            dbSendPreparedQuery(db, "INSERT INTO pedigree VALUES (?,?,?)", as.data.frame(ped))
        }
        if (start.transaction) dbCommit(db)
        invisible(samples)
    }, error=function(e) {
        if (start.transaction) dbRollback(db)
        stop(e)
    })
}

#' Get a list of all samples.
#'
#' @param db database connection
#' @param include.pedigree whether to return pedigree information in the
#' sample table.
#' @param where list of key/value pairs to specify conditions for sample
#' table columns.
#'
#' @return table of sample information
get.all.samples <- function(db, where=NULL, include.pedigree=TRUE) {
    sql <- "SELECT * FROM Sample"
	if (include.pedigree) {
        sql <- paste(sql, "LEFT JOIN Pedigree ON Sample.id = Pedigree.sample_id")
    }
	if (!is.null(where)) {
		sql <- paste(sql, "WHERE", paste(where, collapse=" AND "))
	}
    samples <- dbGetQuery(db, sql)
    invisible(.convert.samples(samples))
}

#' Get a list of the samples in the set with the specified ID, or by a
#' vector of sample IDs.
#'
#' @param ids The sample set ID, or a list of sample IDs.
#' @param db the datbase connection.
#' @param include.pedigree whether to return pedigree information in the
#' sample table.
#' @param where list of key/value pairs to specify conditions for sample
#' table columns.
#'
#' @return data frame with sample information.
get.samples <- function(ids, db, where=NULL, include.pedigree=TRUE) {
    if (length(ids) == 1) {
		if (include.pedigree) {
	        sql <- .wrap("SELECT Sample.*, mother_id, father_id FROM Sample JOIN SampleXSampleSet ON
	            Sample.id = SampleXSampleSet.sample_id LEFT JOIN Pedigree ON Sample.id =
	            Pedigree.sample_id WHERE SampleXSampleSet.set_id = ?")
	    }
	    else {
	        sql <- .wrap("SELECT Sample.* FROM Sample JOIN SampleXSampleSet ON Sample.id =
	            SampleXSampleSet.sample_id WHERE SampleXSampleSet.set_id = ?")
	    }
	}
	else {
		if (include.pedigree) {
			sql <- .wrap("SELECT Sample.*, mother_id, father_id FROM Sample LEFT JOIN Pedigree ON 
				Sample.id = Pedigree.sample_id WHERE Sample.id = ?")
		}
		else {
			sql <- "SELECT * from Sample WHERE id = ?"
		}
	}
    samples <- dbGetPreparedQuery(db, sql, data.frame(ids))
    invisible(.convert.samples(samples))
}

# RSQLite doesn't handle automatic conversion of boolean or date fields
.convert.samples <- function(s) {
	for (col in c("cell_line", "background", "taxon", "source", "description")) {
		w <- !is.na(s[,col]) & s[,col]==""
		if (any(w)) {
			s[w,col] <- NA
		}
	}
    s$impute <- as.logical(s$impute)
    w <- is.na(s$impute)
    if (any(w)) {
        s[w, "impute"] <- FALSE
    }
	# TODO: convert time_added to POSIX datetime
    s
}

#' Import genotype data.
#'
#' \code{genotypes} must be a data.frame with at least three columns. The first two columns specify
#' the sample and marker. Each may be either numeric, in which case it is interpreted as an ID, or
#' character, in which case it is interpreted as a name. The third column is the genotype call.
#' Subsequent columns are interpreted as platform-specific data.
#'
#' Genotypes are .recoded into a bit vector (stored as an integer), where each bit represents one of
#' the possible alleles. Missing information is coded as 0.
#'
#' @param genotypes Genotype data. Either a path to a tab-delimited file (supports gzipped files)
#' or a data.frame. See Details.
#' @param db The database connection.
#' @param N.call The genotype used for missing information (no-call). Any other genotype that is
#' not the two possible alleles or the N.call value is considered a heterozygous call.
#' @param platform The name of the table in which to store platform-specific data, if applicable.
#' @param marker.set.id ID of a marker set. If NULL, marker IDs will be selected by name.
#' @param sample.set.id ID of a sample set. If NULL, sample IDs will be selected by name.
import.genotypes <- function(genotypes, db, N.call="N", platform=NULL, marker.set.id=NULL, sample.set.id=NULL) {
	if (is.character(genotypes)) {
		if (file_ext(genotypes) == "gz") {
			genotypes <- gzfile(genotypes)
		}
		genotypes <- read.table(genotypes, sep="\t", header=T, stringsAsFactors=F, comment.char="")
	}
	
	# Convert sample names to IDs, if necessary
    if (is.character(genotypes[,1])) {
        if (is.null(sample.set.id)) {
            sample.names <- unique(genotypes[,1])
            samples <- dbGetPreparedQuery(db, "SELECT name, id FROM Sample WHERE name = ?",
                data.frame(sample.names))
        }
		else {
			samples <- get.samples(sample.set.id, db)[,c("name","id")]
		}
		m <- match(genotypes[,1], samples[,1])
		if (any(is.na(m))) {
			genotypes <- genotypes[!is.na(m),]
		}
		genotypes[,1] <- samples[m[!is.na(m)],2]
    }

	if (is.character(genotypes[,2])) {
		if (is.null(marker.set.id)) {
			marker.names <- unique(genotypes[,2])
			markers <- dbGetPreparedQuery(db, 
				.wrap("SELECT name, id, alleles from MarkerXMarkerSet 
					join Marker on MarkerXMarkerSet.marker_id = Marker.id 
					WHERE MarkerXMarkerSet.name = ?"),
				data.frame(marker.names))
		}
		else {
			markers <- get.markers(marker.set.id, db)[,c("name","id","alleles")]
		}
        m <- match(genotypes[,2], markers[,1])
		if (any(is.na(m))) {
			genotypes <- genotypes[!is.na(m),]
		}
        genotypes[,2] <- markers[m[!is.na(m)],2]
	}
	
	# .recode genotypes
    recoded <- rep(-1, nrow(genotypes))
    for (a in unique(markers[,3])) {
        w <- genotypes[,2] %in% markers[markers[,3] == a, 2]
        recoded[w] <- .recode(genotypes[w,3], strsplit(a, '')[[1]], N.call)
    }

    tryCatch({
        dbBegin(db)
        dbSendPreparedQuery(db, "INSERT INTO Genotype VALUES (?,?,?)",
            as.data.frame(cbind(genotypes[,1:2], recoded)))

        if (!is.na(platform) && ncol(genotypes) > 3) {
            if (dbExistsTable(db, platform)) {
                dbSendPreparedQuery(db, paste("INSERT INTO", platform, "VALUES (", paste(rep("?",
                    ncol(genotypes) - 1), collapse=","), ")"), genotypes[,-3])
            }
            else {
                stop(paste("Platform-specific data table does not exist", platform))
            }
        }
        dbCommit(db)
    }, error=function(e) {
        dbRollback(db)
        stop(e)
    })
}

#' Fetch genotype data for arbitrary sets of markers and samples.
#'
#' @param markers A marker set ID or a vector of marker IDs.
#' @param db The database connection.
#' @param samples A sample set ID, a list of sample IDs, or NULL for all samples.
#' @param as.matrix Whether the result should be converted from a three-column data frame
#' (marker ID, sample ID, call) to a MxN matrix.
#'
#' @return
#' A MxN data frame of (recoded) genotypes, where M = markers and N = samples. Row names
#' are marker IDs and colnames are sample IDs.
get.genotypes <- function(markers, db, samples=NULL, as.matrix=TRUE) {
    tryCatch({
        dbBegin(db)

        sql <- "SELECT Genotype.marker_id, Genotype.sample_id, call FROM Genotype"
        where <- NULL
        drop <- NULL

        if (length(markers) > 1) {
            dbSendQuery(db, "CREATE TEMPORARY TABLE TempMarkers (id PRIMARY KEY)")
            dbSendPreparedQuery(db, "INSERT INTO TempMarkers VALUES (?)", data.frame(markers))
            sql <- paste(sql, "JOIN TempMarkers ON Genotype.marker_id = TempMarkers.id")
            drop <- c("DROP TABLE TempMarkers")
        }
        else if (!is.null(markers)) {
            sql <- paste(sql, "JOIN MarkerXMarkerSet ON Genotype.marker_id = MarkerXMarkerSet.marker_id")
            where <- c(where, paste("MarkerXMarkerSet.set_id =", markers))
        }

        if (length(samples) > 1) {
            dbSendQuery(db, "CREATE TEMPORARY TABLE TempSamples (id PRIMARY KEY)")
            dbSendPreparedQuery(db, "INSERT INTO TempSamples VALUES (?)", data.frame(samples))
            sql <- paste(sql, "JOIN TempSamples ON Genotype.sample_id = TempSamples.id")
            drop <- c("DROP TABLE TempSamples")
        }
        else if (is.list(samples)) {
            where <- c(where, paste("Genotype.sample_id =", samples[[1]]))
        }
        else if (!is.null(samples)) {
            sql <- paste(sql, "JOIN SampleXSampleSet ON Genotype.sample_id = SampleXSampleSet.sample_id")
            where <- c(where, paste("SampleXSampleSet.set_id =", samples))
        }

        if (!is.null(where)) {
            sql <- paste(sql, "WHERE", paste(where, collapse=" AND "))
        }

        geno <- dbGetQuery(db, sql)

        if (!is.null(drop)) {
            for (d in drop) {
                dbSendQuery(db, d)
            }
        }
    }, error=function(e) {
        stop(e)
    }, finally={
        dbRollback(db)
    })
    if (as.matrix) {
        geno <- .my.dcast(geno, marker_id~sample_id, value.var="call")
    }
    invisible(geno)
}

#' Fetch platform-specific data.
#'
#' @param marker.set.id The assay ID.
#' @param platform The name of the platform-specific data table.
#' @param db The database connection.
#' @param samples The set of samples, or all samples if NULL.
#' @param as.matrix Whether the result should be converted from a four column data frame to a list
#' of two MxN matricies.
#' @param format If as.matrix=TRUE, controls whether rows are samples and columns are markers
#' (format=1) or vice-versa (format=2).
#'
#' @return
#' A list of data frames. Each is a MxN data frame of a platform-specific variable,
#' where M = markers and N = samples.
get.platform.data <- function(marker.set.id, platform, db, samples=NULL, as.matrix=TRUE, format=1) {
    sql <- .wrap(paste("SELECT ", platform, ".* FROM ", platform, " JOIN MarkerXMarkerSet ON ", 
        platform, ".marker_id = MarkerXMarkerSet.marker_id", sep=""))
    if (is.null(samples)) {
        sql <- paste(sql, "WHERE")
        params <- marker.set.id
    }
    else if (length(samples) == 1) {
        sql <- paste(sql, paste("JOIN SampleXSampleSet ON ", platform, 
            ".sample_id = SampleXSampleSet.sample_id WHERE SampleXSampleSet.set_id = ? AND ", sep=""))
        params <- data.frame(samples, marker.set.id)
    }
    else {
        sql <- paste(sql, paste("WHERE ", platform, ".sample_id = ? AND ", sep=""))
        params <- data.frame(samples, rep(marker.set.id, length(samples)))
    }
    sql <- paste(sql, "MarkerXMarkerSet.set_id = ?")
    data <- dbGetPreparedQuery(db, sql, params)
    if (as.matrix) {
        n <- colnames(data)[3:ncol(data)]
        if (format == 1) {
            data <- lapply(3:ncol(data), function(i) {
                .my.dcast(data, sample_id~marker_id, value.var=colnames(data)[i])
            })
        }
        else {
            data <- lapply(3:ncol(data), function(i) {
                .my.dcast(data, marker_id~sample_id, value.var=colnames(data)[i])
            })
        }
        names(data) <- n
    }
    invisible(data)
}

#' Insert assay into database.
#'
#' @param assay CLASPAssay object.
#' @param db database connection.
#'
#' @return the assay updated with an ID column.
insert.assay <- function(assay, db) {
	now <- Sys.time()
	dbSendPreparedQuery(db, .wrap("INSERT INTO assay (name, marker_set_id, sample_set_id, ob_name, 
		ob_marker_set_id, ob_sample_set_id, time_added) VALUES (?,?,?,?,?,?,?)"), cbind(assay, now))
	assay.id <- dbGetQuery(db, "SELECT last_insert_rowid() FROM assay")[1,1]
	assay <- cbind(id=assay.id, assay)
	invisible(assay)
}

#' Fetch assay from the database.
#'
#' @param id assay ID.
#' @param db database connection.
#' @param load logical; if true, sample, marker and genotype data are loaded from the database.
#'
#' @return a list containing fields from the Assay table. If load=TRUE, the list also contains
#' loaded sample, marker and genotype data AND the object is designated as a CLASPAssay object.
get.assay <- function(id, db, load=TRUE) {
	assay <- as.list(dbGetQuery(db, "SELECT * FROM assay WHERE id = ?", id)[1,])
	if (load) {
		geno <- get.genotypes(assay$marker_set_id, db, assay$sample_set_id)
		samples <- get.samples(assay$sample_set_id, db)
		ob.geno <- get.genotypes(assay$marker_set_id, db, assay$sample_set_id)
		ob.samples <- get.samples(assay$sample_set_id, db)
		assay$inbred <- list(
			markers=get.marker.set(assay$marker_set_id, db), samples=samples, 
			geno=geno, distmat=.distance.matrix(.pairwise.distances(geno)), samples$name)
		assay$outbred <- list(
			markers=get.marker.set(assay$marker_set_id, db),
			samples=ob.samples, geno=ob.geno, distmat=.orthoganal.distance.matrix(
				ob.geno, geno[match(rownames(ob.geno), rownames(geno)),], 
				dnames=list(ob.samples$name, samples$name))
		)
		class(assay) <- "CLASPAssay"
	}
	invisible(assay)
}
