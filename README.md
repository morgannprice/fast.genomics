fast.genomics is a fast comparative genome browser for thousands of
representative genomes (1 per genus) of bacteria and archeaa.
fast.genomics also has "sub databases" for each order of bacteria and
archaea. The sub databases represent ever species that has a high
quality genome, and include up to 10 genomes per species. This
supports fast browsing for every group of prokaryotes.

See
http://fast.genomics.lbl.gov/
or the [preprint](https://doi.org/10.1101/2023.08.23.554478)

# System Requirements

These scripts should work on any Linux system, and would probably work
on other Unix or MacOS as well. All of the code is written in perl (I
am using perl v5.16.3.). The web site uses the common gateway
interface (CGI).

Some of the scripts require at least 64 GB of RAM. In particular, the
main database is designed for the mmseqs2 index to be stored in memory
(53 GB for 4,788 representative genomes).

# Downloads

You can download the current database as a protein fasta file and a
SQLite3 database, see bottom of

http://fast.genomics.lbl.gov/

# Installation

After checking out the code, see the instructions in SETUP.  Also see
bin/fetchAssemblies.pl and bin/build.pl, which together build the
database by downloading information from RefSeq and Genbank. Once you
have the main database built, you can use bin/buildSubDbs.pl to build
a "sub db" for each order.
