#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.2';
my $name = 'run_phylogenies.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Quick phylogenetic pipeline based on MAFFT, TRIMAL and IQ-TREE
REQUIREMENTS	MAFFT - https://mafft.cbrc.jp/alignment/software/ (doi: 10.1093/bioinformatics/bty121)
		TRIMAL - http://trimal.cgenomics.org/ (doi: 10.1093/bioinformatics/btp348)
		IQTREE - http://www.iqtree.org/ (doi: 10.1093/molbev/msu300)
		
USAGE		run_phylogenies.pl -t 10 -f *.fasta

OPTIONS:
-t (--threads)	Number of threads to use [Default: 10]
-f (--fasta)	Files in Multifasta format

## MAFFT, TRIMAL and IQ-TREE command lines can be edited at the bottom of this script
OPTIONS
die "$options\n" unless @ARGV;

my $threads = 10;
my @files;
GetOptions(
	't|threads=i'=>\$threads,
	'f|fasta=s@{1,}'=>\@files
);

## Checking if programs are found
my $mafft = `command -v mafft`; chomp $mafft; if ($mafft eq ''){print "\nERROR: Cannot find MAFFT. Please install MAFFT in your \$PATH\n\n"; exit;}
my $trimal = `command -v trimal`; chomp $trimal; if ($trimal eq ''){print "\nERROR: Cannot find trimAl. Please install trimAl in your \$PATH\n\n"; exit;}
my $iqtree = `command -v iqtree`; chomp $iqtree; if ($iqtree eq ''){print "\nERROR: Cannot find IQ-TREE. Please install IQ-TREE in your \$PATH\n\n"; exit;}

## Running phylogenies
while (my $fasta = shift@files){
	my ($ext) = $fasta =~ /.(\w+)$/;
	$fasta =~ s/.\w+$//;
	system "echo Aligning $fasta with MAFTT...";
	system "mafft --thread $threads $fasta.$ext > $fasta.aln";
	system "echo Trimming $fasta alignment with TrimAL...";
	system "trimal -in $fasta.aln -out ${fasta}_trimmed.tfa -htmlout $fasta.html -automated1";
	system "Computing $fasta phylogenetic tree with IQtree...";
	system "iqtree -s ${fasta}_trimmed.tfa -nt $threads -bb 1000";
}
