CREATE TABLE Genome (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    taxon       TEXT,
    version     TEXT
);
CREATE TABLE Array (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT,
    platform        TEXT
);
CREATE TABLE ArrayXGenome (
    array_id        INTEGER,
    genome_id       INTEGER,
    marker_set_id   INTEGER,
    FOREIGN KEY(array_id) REFERENCES Array(id),
    FOREIGN KEY(genome_id) REFERENCES Genome(id),
    FOREIGN KEY(marker_set_id) REFERENCES MarkerSet(id)
);
CREATE TABLE MarkerSet (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT,
    description TEXT
);
CREATE TABLE Marker (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    genome_id   INTEGER,
    chromosome  INTEGER,
    position_bp INTEGER,
    position_cm REAL,
    alleles     TEXT,
    FOREIGN KEY(genome_id) REFERENCES Genome(id)
);
CREATE TABLE MarkerXMarkerSet (
    marker_id   INTEGER,
    set_id      INTEGER,
    name        TEXT,
    FOREIGN KEY(marker_id) REFERENCES Marker(id),
    FOREIGN KEY(set_id) REFERENCES MarkerSet(id)
);
CREATE TABLE SampleSet (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT,
    description TEXT
);
CREATE TABLE Sample (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT,
    cell_line   TEXT,
    background  TEXT,
    type        TEXT,
    taxon       TEXT,
    source      TEXT,
    description TEXT,
	impute		BOOLEAN,
    array_id    INTEGER,
    time_added  DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY(array_id) REFERENCES Array(id)
);
CREATE TABLE SampleXSampleSet (
    sample_id   INTEGER,
    set_id      INTEGER,
    FOREIGN KEY(set_id) REFERENCES SampleSet(id),
    FOREIGN KEY(sample_id) REFERENCES Sample(id)
);
CREATE TABLE Pedigree (
    sample_id   INTEGER,
    mother_id   INTEGER,
    father_id   INTEGER,
    FOREIGN KEY(sample_id) REFERENCES Sample(id),
    FOREIGN KEY(mother_id) REFERENCES Sample(id),
    FOREIGN KEY(father_id) REFERENCES Sample(id)
);
CREATE TABLE Genotype (
    sample_id   INTEGER,
    marker_id   INTEGER,
    call        INTEGER,
    FOREIGN KEY(marker_id) REFERENCES Marker(id),
    FOREIGN KEY(sample_id) REFERENCES Sample(id)
);
CREATE TABLE Illumina (
    sample_id   INTEGER,
    marker_id   INTEGER,
    X           REAL,
    Y           REAL,
    FOREIGN KEY(marker_id) REFERENCES Marker(id),
    FOREIGN KEY(sample_id) REFERENCES Sample(id)
);
CREATE TABLE Assay (
	id                  INTEGER PRIMARY KEY AUTOINCREMENT,
	name                TEXT,
	marker_set_id       INTEGER, 
	sample_set_id       INTEGER, 
	ob_name             TEXT, 
	ob_marker_set_id    INTEGER,
	ob_sample_set_id    INTEGER, 
	time_added          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY(marker_set_id) REFERENCES MarkerSet(id),
    FOREIGN KEY(sample_set_id) REFERENCES SampleSet(id),
    FOREIGN KEY(ob_marker_set_id) REFERENCES MarkerSet(id),
    FOREIGN KEY(ob_sample_set_id) REFERENCES SampleSet(id)
);
CREATE UNIQUE INDEX Genome_taxon_version ON Genome (taxon, version);
CREATE UNIQUE INDEX Marker_genome_location ON Marker (genome_id, chromosome, position_bp);
CREATE INDEX Marker_location ON Marker (chromosome, position_bp);
CREATE INDEX MarkerXMarkerSet_marker_set ON MarkerXMarkerSet(marker_id, set_id);
CREATE INDEX MarkerXMarkerSet_set_marker ON MarkerXMarkerSet(set_id, marker_id);
CREATE INDEX SampleXSampleSet_sample_set ON SampleXSampleSet(sample_id, set_id);
CREATE INDEX SampleXSampleSet_set_sample ON SampleXSampleSet(set_id, sample_id);
CREATE INDEX Genotype_marker_sample ON Genotype (marker_id, sample_id);
CREATE INDEX Genotype_sample_marker ON Genotype (sample_id, marker_id);
CREATE INDEX Illumina_marker_sample ON Illumina (marker_id, sample_id);
CREATE INDEX Illumina_sample_marker ON Illumina (sample_id, marker_id);
