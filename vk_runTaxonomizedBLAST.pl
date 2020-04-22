#!/usr/bin/perl
## Pombert Lab, IIT 2016

## Requirements:
## BLAST 2.2.28+ or later
## NCBI taxonomony database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)
## NCBI NR/NT databases (ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
## NOTE: the NCBI taxdb, nr and nt databases can be downloaded with the update_blastdb.pl from NCBI (http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl)
## The BLASTDB variable must be set in the environmental variables: export BLASTDB=/path/to/NCBI/TaxDB:/path/to/NCBI/NR:/path/to/NCBI/NT

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "
USAGE = perl runTaxonomizedBLAST.pl [options]

OPTIONS:
--type		blastn, blastx, blastp... [default = blastp]
--db		nt, nr or custom subset [default = nr]
--threads	2, 4, 8 or more [default = 2]
--evalue	1e-05, 1e-10 or other [default = 1e-05]
--culling	culling limit [default = 1]
--query		fasta file to be queried
";

die "$usage\n" unless@ARGV;

## Defining options
my $blast_type = 'blastp';
my $db = 'nr';
my $threads = '2';
my $evalue = '1e-05';
my $culling = '1';
my $query;

GetOptions(
    'type=s' => \$blast_type,
	'db=s' => \$db,
	'threads=i' => \$threads,
	'evalue=s' => \$evalue,
	'culling=i' => \$culling,
	'query=s' => \$query,
);

## Running BLAST
system "$blast_type -num_threads $threads -query $query -db $db -evalue $evalue -culling_limit $culling -outfmt '6 qseqid sseqid qstart qend pident length bitscore evalue staxids sscinames sskingdoms sblastnames' -out $query.$blast_type";

