CREATE TABLE KmerFirst (
  kmer TEXT NOT NULL,
  seek BIGINT NOT NULL,
  PRIMARY KEY(kmer,seek)
);

CREATE TABLE KmerLast (
  kmer TEXT NOT NULL,
  seek BIGINT NOT NULL,
  PRIMARY KEY(kmer,seek)
);
