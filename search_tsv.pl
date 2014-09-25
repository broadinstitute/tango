# Given read names, specified in a tsv file as the 1st field, obtain actual reads 
# (optional to output as FASTA or FASTQ format) by parsing through a given set of 
# fq files separated by comma in the commandline

#!/usr/bin/perl
use strict;
use Getopt::Long;

my $query = shift;
my $target = shift;
my $num_reads = 0;

# ------------- read through query file -------------
my %rname_hash;
my %rstr_hash;
open (QUERY, $query) or die $!;
while (my $line = <QUERY>) {
	next if ($line =~ /^\s*$/); #skip blank lines
	chomp($line);
	my @entries = split ('\t', $line);
	$rname_hash{$entries[0]} = 1;	
	$rstr_hash{$entries[1]} = 1;
	++ $num_reads;
}
close (QUERY);
print "\n\tnum_reads in $query: $num_reads\n\n";
$num_reads = 0;

# ---------- search in target file ------------------
 
open (TARGET, $target) or die $!;
while (my $line = <TARGET>) {
	++ $num_reads;
	if ($num_reads % 1000000 == 0) {
#		print $num_reads."\n";
	}
	next if ($line =~ /^\s*$/); #skip blank lines
	chomp($line);
	my @entries = split ('\t', $line);
	my $header = $entries[0];
	my $read = $entries[1];
	if (exists $rstr_hash{$read}) {
		if (exists $rname_hash{$header}) {
			print "\nTARGET\n";
		}
		print "@entries\n";	
	}	
}
close (TARGET);
print "\n\tnum_reads in $target: $num_reads\n\n";


################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Given read names specified in a tsv file as the 1st field and a list of fastq files,\n";
		print "\tsearching the reads (optional to output as FASTA or FASTQ format) in the fastq files\n";
		print "\tand generate the read statistics\n\n";
		print "usage: ./$0 {-ihit [hits.tsv] -ifqlist [1.fq],[2.fq],[3.fq] \n\n";
		print "-h: help; for cmd options\n";
		print "-silent: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-ihit: input hit file in tsv format with 1st field being read name\n";
		print "-ifqlist: input fastq file list\n";
		print "-o: opt; when specified, the reads will be printed out to this file (.fa .fasta .fq .fastq)\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

 