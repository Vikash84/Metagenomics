#!/usr/bin/perl
## Pombert Lab, 2018
## This script generates a distribution of sequences per species, genus, family and so forth from taxonomised BLAST output files.
## This script was created to handle megablast analyses of nanopore 16S amplicon sequencing.
## Because of the error-rate of the nanopore sequencing, identification at the species/subspecies level can be ambiguous and should not be taken at face value.
## Values at the genus level should be more accurate given the observed error-rate of the technology.
## Requires ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
## v 0.1

use strict;
#use warnings;
use Getopt::Long qw(GetOptions);

### Defining options
my $options = <<'OPTIONS';

EXAMPLE: taxid_dist.pl -n TaxDumps/nodes.dmp -a TaxDumps/names.dmp -b Examples/*.blastn -e 1e-75 -h 1

NOTE: requires the following BLAST format: -outfmt '6 qseqid sseqid pident length bitscore evalue staxids sskingdoms sscinames sblastnames'

-n (--nodes)	NCBI nodes.dmp file 
-a (--names)	NCBI names.dmp
-b (--blast)	NCBI blast output file(s) in oufmt 6 format
-e (--evalue)	evalue cutoff [Default: 1e-75]
-h (--hits)	Number of BLAST hits to keep; top N hits [Default: 1]

OPTIONS
die "$options" unless @ARGV;

my $node;
my $name;
my @blast = ();
my $evalue = 1e-75;
my $hits = 1;

GetOptions(
	'n|nodes=s' => \$node,
	'a|names=s' => \$name,
	'b|blast=s@{1,}' => \@blast,
	'e|value=s' => \$evalue,
	'h|hits=i' => \$hits
);

## Initializing taxids -> names database
my %taxid;
open NAMES, "<$name";
system "echo Initializing taxonomic IDs...";
while (my $line = <NAMES>){
	chomp $line; $line =~ s/\t\|//g;
	if ($line =~ /scientific name/){
		my @columns = split("\t", $line);
		$taxid{$columns[0]} = $columns[1];
	}
}

### Initializing taxonomic databases
my %phylum; my %norank; my %class; my %order; my %family; my %genus; my %species; my %subspecies;
open NODES, "<$node";
system "echo Initializing taxonomic databases...";
while (my $line = <NODES>){
	chomp $line; $line =~ s/\t\|//g;
	my @columns = split("\t", $line);	## $columns[0] = taxid; $columns[1] = parent tax_id; $columns[2] = rank
	if ($columns[2] eq "no rank"){$norank{$columns[0]}[0] = $taxid{$columns[0]}; $norank{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "subspecies"){$subspecies{$columns[0]}[0] = $taxid{$columns[0]}; $subspecies{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "species"){$species{$columns[0]}[0] = $taxid{$columns[0]}; $species{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "genus"){$genus{$columns[0]}[0] = $taxid{$columns[0]}; $genus{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "family"){$family{$columns[0]}[0] = $taxid{$columns[0]}; $family{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "order"){$order{$columns[0]}[0] = $taxid{$columns[0]}; $order{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "class"){$class{$columns[0]}[0] = $taxid{$columns[0]};	$class{$columns[0]}[1] = $columns[1];}
	elsif ($columns[2] eq "phylum"){$phylum{$columns[0]}[0] = $taxid{$columns[0]}; $phylum{$columns[0]}[1] = $columns[1];}
}

## Working on BLAST taxonomized outfmt6 files; treating each file as independent datasets
## -outfmt '6 qseqid sseqid pident length bitscore evalue staxids sskingdoms sscinames sblastnames'
my %bhits; my %bphylum; my %bnorank; my %bclass; my %border; my %bfamily; my %bgenus; my %bspecies; my %bsubspecies;
my @col; my $sub; my $spp; my $gen; my $fam; my $ord; my $class; my $nrk; my $phy;
while (my $blast = shift@blast){
	system "echo Parsing $blast...";
	open BLAST, "<$blast";
	open SUB, ">$blast.subspecies";	open SPP, ">$blast.species"; open GEN, ">$blast.genus";	open FAM, ">$blast.family";
	open ORD, ">$blast.order"; open CLAS, ">$blast.class"; open PHY, ">$blast.phylum"; open NORK, ">$blast.norank";
	%bhits = (); %bphylum = (); %bnorank = (); %bclass = (); %border = (); %bfamily = ();	%bgenus = (); %bspecies = (); %bsubspecies = ();
	$sub = 0; $spp = 0; $gen = 0; $fam = 0; $ord = 0; $class = 0; $nrk = 0; $phy = 0;
	while (my $line = <BLAST>){
		chomp $line;
		@col = split("\t", $line);	## $col[0] = qseqid; $col[1] = sseqid; $col[2] = pident; $col[3] = length; $col[4] = bitscore; 
								## $col[5] = evalue; $col[6] = staxids; $col[7] = sskingdoms; $col[8] = sscinames; $col[9] = sblastnames
		if ($col[6] =~ /^(\d+);/){$col[6] = $1;} ## Searching for multiple taxids, keeping only the 1st one
		if ($col[5] <= $evalue){
			if ((exists $bhits{$col[0]}) && ($bhits{$col[0]} >= $hits)){next;}
			elsif ((exists $bhits{$col[0]}) && ($bhits{$col[0]} < $hits)){
				$bhits{$col[0]}++;
				if (exists $norank{$col[6]}){$bnorank{$col[6]} += 1; $nrk++;}
				elsif (exists $subspecies{$col[6]}){subspecies();}
				elsif (exists $species{$col[6]}){species();}
				else {print "ERROR: $col[6] taxonomic level is not a species, subspecies, or incertae sedis (no rank)\n";}
			}
			else{
				$bhits{$col[0]} = 1;
				if (exists $norank{$col[6]}){$bnorank{$col[6]} += 1; $nrk++;}
				elsif (exists $subspecies{$col[6]}){subspecies();}
				elsif (exists $species{$col[6]}){species();}
				else {print "ERROR: $col[6] taxonomic level is not a species, subspecies, or incertae sedis (no rank)\n";}
			}
		}
	}
	my $size;
	print NORK "No rank\tTaxIF\tNumber\tPercent (total = $nrk)\n"; foreach (sort {$bnorank{$b} <=> $bnorank{$a}} keys %bnorank){my $av = sprintf("%.2f%%", ($bnorank{$_}/$nrk)*100); print NORK "$taxid{$_}\t$_\t$bnorank{$_}\t$av\n";}
	print SUB "Subspecies\tTaxIF\tNumber\tPercent (total = $sub)\n"; foreach (sort {$bsubspecies{$b} <=> $bsubspecies{$a}} keys %bsubspecies){my $av = sprintf("%.2f%%", ($bsubspecies{$_}/$sub)*100); print SUB "$taxid{$_}\t$_\t$bsubspecies{$_}\t$av\n";}
	print SPP "Species\tTaxIF\tNumber\tPercent (total = $spp)\n"; foreach (sort {$bspecies{$b} <=> $bspecies{$a}} keys %bspecies){my $av = sprintf("%.2f%%", ($bspecies{$_}/$spp)*100); print SPP "$taxid{$_}\t$_\t$bspecies{$_}\t$av\n";}
	print GEN "Genus\tTaxIF\tNumber\tPercent (total = $gen)\n"; foreach (sort {$bgenus{$b} <=> $bgenus{$a}} keys %bgenus){my $av = sprintf("%.2f%%", ($bgenus{$_}/$gen)*100); print GEN "$taxid{$_}\t$_\t$bgenus{$_}\t$av\n";}
	print FAM "Family\tTaxIF\tNumber\tPercent (total = $fam)\n"; foreach (sort {$bfamily{$b} <=> $bfamily{$a}} keys %bfamily){my $av = sprintf("%.2f%%", ($bfamily{$_}/$fam)*100); print FAM "$taxid{$_}\t$_\t$bfamily{$_}\t$av\n";}
	print ORD "Order\tTaxIF\tNumber\tPercent (total = $ord)\n"; foreach (sort {$border{$b} <=> $border{$a}} keys %border){my $av = sprintf("%.2f%%", ($border{$_}/$ord)*100); print ORD "$taxid{$_}\t$_\t$border{$_}\t$av\n";}
	print CLAS "Class\tTaxIF\tNumber\tPercent (total = $class)\n"; foreach (sort {$bclass{$b} <=> $bclass{$a}} keys %bclass){my $av = sprintf("%.2f%%", ($bclass{$_}/$class)*100); print CLAS "$taxid{$_}\t$_\t$bclass{$_}\t$av\n";}
	print PHY "Class\tTaxIF\tNumber\tPercent (total = $phy)\n"; foreach (sort {$bphylum{$b} <=> $bphylum{$a}} keys %bphylum){my $av = sprintf("%.2f%%", ($bphylum{$_}/$phy)*100); print PHY "$taxid{$_}\t$_\t$bphylum{$_}\t$av\n";}
}
## Subroutines
sub subspecies{
	$bsubspecies{$col[6]} += 1; $sub++;
	$bspecies{$subspecies{$col[6]}[1]} += 1; $spp++;
	$bgenus{$species{$subspecies{$col[6]}[1]}[1]} += 1; $gen++;
	$bfamily{$genus{$species{$subspecies{$col[6]}[1]}[1]}[1]} += 1; $fam++;
	$border{$family{$genus{$species{$subspecies{$col[6]}[1]}[1]}[1]}[1]} += 1; $ord++;
	$bclass{$order{$family{$genus{$species{$subspecies{$col[6]}[1]}[1]}[1]}[1]}[1]} += 1; $class++;
	$bphylum{$class{$order{$family{$genus{$species{$subspecies{$col[6]}[1]}[1]}[1]}[1]}[1]}[1]} += 1; $phy++;
}
sub species{
	$bspecies{$col[6]} += 1; $spp++;
	$bgenus{$species{$col[6]}[1]} += 1; $gen++;
	$bfamily{$genus{$species{$col[6]}[1]}[1]} += 1; $fam++;
	$border{$family{$genus{$species{$col[6]}[1]}[1]}[1]} += 1; $ord++;
	$bclass{$order{$family{$genus{$species{$col[6]}[1]}[1]}[1]}[1]} += 1; $class++;
	$bphylum{$class{$order{$family{$genus{$species{$col[6]}[1]}[1]}[1]}[1]}[1]} += 1; $phy++;
}