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

CREATE TABLE Taxon (
  taxon TEXT NOT NULL,
  level TEXT NOT NULL,
  parent TEXT NOT NULL,
  nGenomes INT NOT NULL,
  PRIMARY KEY (taxon,level)
);
CREATE INDEX TaxonLevel ON Taxon (level, taxon);
CREATE INDEX TaxonParent ON Taxon (parent);
