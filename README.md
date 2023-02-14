fast.genomics is a fast comparative genome browser for thousands of
representative genomes (1 per genus) of bacteria and archeaa. See

http://fast.genomics.lbl.gov/

# System Requirements

These scripts should work on any Linux system, and would probably work
on other Unix or MacOS as well. All of the code is written in perl (I
am using perl v5.16.3.). The web site uses the common gateway
interface (CGI).

# Downloads

You can download the current database as a protein fasta file and a
SQLite3 database, see bottom of

http://fast.genomics.lbl.gov/

# Installation

After checking out the code, see the instructions in SETUP.  Also see
bin/fetchAssemblies.pl and bin/build.pl, which together build the
database by downloading information from RefSeq and Genbank.
