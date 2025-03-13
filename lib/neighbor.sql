/* This is the schema for both the top-level database, which stores
   the representative genomes for each genus, and for each "subdb",
   which has many more genomes, but is restricted to a taxon (usually
   one subdb per order).

   A few fields or tables are used only in the top-level database
   or only in the subdbs.
*/

CREATE TABLE Genome (
  gid TEXT PRIMARY KEY, /* assemblyId in the downloaded tables */
  gtdbDomain TEXT NOT NULL,
  gtdbPhylum TEXT NOT NULL,
  gtdbClass TEXT NOT NULL,
  gtdbOrder TEXT NOT NULL,
  gtdbFamily TEXT NOT NULL,
  gtdbGenus TEXT NOT NULL,
  gtdbSpecies TEXT NOT NULL,
  strain TEXT NOT NULL,
  gtdbAccession TEXT NOT NULL,
  assemblyName TEXT NOT NULL,
  ncbiTaxonomy TEXT NOT NULL,
  nGenes INT NOT NULL,
  nProteins INT NOT NULL /* #protein-coding genes */
);
CREATE INDEX 'GenomeDomain' ON Genome (gtdbDomain);
CREATE INDEX 'GenomePhylum' ON Genome (gtdbPhylum);
CREATE INDEX 'GenomeClass' ON Genome (gtdbClass);
CREATE INDEX 'GenomeOrder' ON Genome (gtdbOrder);
CREATE INDEX 'GenomeFamily' ON Genome (gtdbFamily);
CREATE INDEX 'GenomeGenus' ON Genome (gtdbGenus);
CREATE INDEX 'GenomeSpecies' ON Genome (gtdbSpecies);
CREATE INDEX 'GenomeStrain' ON Genome (strain);
CREATE INDEX 'GenomeAcc' ON Genome (gtdbAccession);
CREATE INDEX 'GenomeAssembly' ON Genome (assemblyName);

/* This table is only built for the top-level database.
   It indicates if a genome is present in the top-level database or not,
   and which SubDB prefix it is in, if any */
CREATE TABLE AllGenome (
  gid TEXT PRIMARY KEY, /* assemblyId in the downloaded tables */
  gtdbDomain TEXT NOT NULL,
  gtdbPhylum TEXT NOT NULL,
  gtdbClass TEXT NOT NULL,
  gtdbOrder TEXT NOT NULL,
  gtdbFamily TEXT NOT NULL,
  gtdbGenus TEXT NOT NULL,
  gtdbSpecies TEXT NOT NULL,
  strain TEXT NOT NULL,
  gtdbAccession TEXT NOT NULL,
  assemblyName TEXT NOT NULL,
  ncbiTaxonomy TEXT NOT NULL,
  nGenes INT NOT NULL,
  nProteins INT NOT NULL, /* #protein-coding genes */
  inTop INT NOT NULL,
  prefix TEXT NOT NULL /* the subdb; empty if not in a subdb */
);
CREATE INDEX 'AllGenomeDomain' ON AllGenome (gtdbDomain);
CREATE INDEX 'AllGenomePhylum' ON AllGenome (gtdbPhylum);
CREATE INDEX 'AllGenomeClass' ON AllGenome (gtdbClass);
CREATE INDEX 'AllGenomeOrder' ON AllGenome (gtdbOrder);
CREATE INDEX 'AllGenomeFamily' ON AllGenome (gtdbFamily);
CREATE INDEX 'AllGenomeGenus' ON AllGenome (gtdbGenus);
CREATE INDEX 'AllGenomeSpecies' ON AllGenome (gtdbSpecies);
CREATE INDEX 'AllGenomeStrain' ON AllGenome (strain);
CREATE INDEX 'AllGenomeAcc' ON AllGenome (gtdbAccession);
CREATE INDEX 'AllGenomeAssembly' ON AllGenome (assemblyName);

CREATE TABLE Scaffold (
  gid TEXT NOT NULL,
  scaffoldId TEXT NOT NULL,
  scaffoldDesc TEXT NOT NULL,
  length INT NOT NULL,
  PRIMARY KEY (gid, scaffoldId)
);

CREATE TABLE Gene (
  gid TEXT NOT NULL,
  scaffoldId TEXT NOT NULL,
  begin INT NOT NULL,
  end INT NOT NULL,
  strand TEXT NOT NULL,
  locusTag TEXT NOT NULL,  /* like AF_RS00005 */
  proteinId TEXT NOT NULL, /* like WP_010877517.1, or the empty string if not protein coding */
  desc TEXT NOT NULL,
  PRIMARY KEY (gid, scaffoldId, begin, end, strand, locusTag)
);
CREATE UNIQUE INDEX 'GeneLocusTag' ON Gene (locusTag);
CREATE INDEX 'GeneProtein' ON Gene (proteinId);

CREATE TABLE Protein (
  proteinId TEXT PRIMARY KEY,
  sequence TEXT NOT NULL
);

/* Cluster is not used in the top-level database */
CREATE TABLE ClusterProtein (
  clusterId TEXT NOT NULL, /* the proteinId of the representative */
  proteinId TEXT NOT NULL,
  PRIMARY KEY (clusterId,proteinId)
);
CREATE INDEX 'ProteinToCluster' ON ClusterProtein (proteinId,clusterId);

/* ClusteringInfo is not used in the top-level database */
CREATE TABLE ClusteringInfo (
  nProteins INT NOT NULL,
  nClusters INT NOT NULL,
  nClusteredAA INT NOT NULL,
  nAA INT NOT NULL
);

/* In the subdbs, only taxa up to the order level are represented */
CREATE TABLE Taxon (
  taxon TEXT NOT NULL,
  level TEXT NOT NULL,
  parent TEXT NOT NULL,
  nGenomes INT NOT NULL,
  PRIMARY KEY (taxon,level)
);
CREATE INDEX TaxonLevel ON Taxon (level, taxon);
CREATE INDEX TaxonParent ON Taxon (parent);

/* SubDbs are listed only in the top-level database */
CREATE TABLE SubDb (
  taxon TEXT NOT NULL PRIMARY KEY,
  level TEXT NOT NULL,
  prefix TEXT NOT NULL,
  nGenomes INT NOT NULL,
  nProteins INT NOT NULL,
  nClusters INT NOT NULL
);

/* All 16S loci (after filtering out short sequences) */
CREATE TABLE All16S (
  gid TEXT NOT NULL,
  locusTag TEXT NOT NULL,
  sequence TEXT NOT NULL,
  PRIMARY KEY (gid,locusTag)
);
CREATE INDEX All16SLocus ON All16S (locusTag, gid);

