#!/usr/bin/perl
## Pombert JF, IIT 2015

## Reorient RNAs in a multifasta file according to their top blast hits.
## RNAs with blast hits are oriented 5' to 3'.
## RNAs with no hits are put into a separate file.
## Require a blastx output file in tabular (outfmt 6) format.

my $usage = "orientRNAs.pl multifasta_file blastx_output desired_file_name";
die "USAGE = $usage\n" unless(@ARGV==3);

open FASTA, "<$ARGV[0]";
open BLAST, "<$ARGV[1]";
open RNA, ">$ARGV[2].oriented.fasta";
open UNK, ">$ARGV[2].unknown.fasta";

## Creating hash of sequences
my %fasta = ();
my @list = ();
my $key = undef;
my $seq = undef;
while (my $line = <FASTA>){
	chomp $line;
	if ((defined$key)&&($line =~ />(.*)$/)){
		my $tmp = $1;
		$fasta{$key}=$seq;
		$key = $tmp;
		push(@list, $key);
		$seq = undef;
	}
	elsif($line =~ />(.*)$/){$key=$1;push(@list, $key);}
	else {$seq .= $line;}
}
$fasta{$key}=$seq;
## End of hash

## Parsing the blast output
my %blast = ();
while (my $line = <BLAST>){
	chomp $line;
	if ($line =~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/){
		my $query = $1;
		my $hit = $2;
		my $id = $3;	## percent identity
		my $len = $4;	## alignment length
		my $mis = $5;	## number of mismatches
		my $gaps = $6; 	## number of gaps
		my $start = $7;	## query start
		my $end = $8;	## query end
		my $tstart = $9;	## target start
		my $tend = $10;	## target end
		my $evalue = $11;
		my $bit	= $12;	## bitscore
		if (exists $blast{$query}){next;}
		else{
			$blast{$query}=$hit;
			if($start <= $end){ ## Looking if RNA is on the forward strand (based on top blast hit)
				print RNA ">$query\n";
				print RNA "$_\n" for unpack '(A60)*', $fasta{$query};
			}
			elsif($start >= $end){ # Looking if RNA is on the reverse strand (based on top blast hit)
				my $rna = reverse($fasta{$query});
				$rna =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
				print RNA ">$query\n";
				print RNA "$_\n" for unpack '(A60)*', $rna;
			}
		}
	}
}

## Writing the sequences without any blast hit into a separate file
while (my $sequence = shift@list){
	chomp $sequence;
	if (exists $blast{$sequence}){next;}
	else{
		print UNK ">$sequence\n";
		print UNK "$_\n" for unpack '(A60)*', $fasta{$sequence};
	}
}
close FASTA; close BLAST; close RNA; close UNK
