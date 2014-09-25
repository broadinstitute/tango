#	@Brief
#		  	Given read names specified in a tsv file as the 1st field and a list of fastq files
#		  	separated by comma in the commandline, searching the reads (optional output as FASTA 
#		  	or FASTQ format if -o) in the fastq files and	generate the read statistics.
#		  	Note: the query read name should not contain the first @ char of fq format.

#!/usr/bin/perl
use strict;
use Getopt::Long;

my %option = (
	h         => '',
	silent 		=> '',
	# I/O setting 
	ihit			=> '', 
	ifqlist		=> '',
	o					=> '',
);

GetOptions(
	"h"				=> \$option{h},
	"silent"	=> \$option{silent},
	"ihit=s"			=> \$option{ihit},
	"ifqlist=s"			=> \$option{ifqlist},
	"o=s"			=> \$option{o},
) || printHelp (); 

if ($option{h}) { printHelp();}

my $ihit = $option{ihit};
my $ifqlist = $option{ifqlist};
my @fqfiles = split (',', $ifqlist);
my $output = $option{o};
my $to_prnt = 0;
my $is_fastq = 0;
if ($output) {
	lc ($output);
	if ($output =~ /\.fq/ || $output =~ /\.fastq/) {
		$is_fastq = 1;
	} else {
		unless ($output =~ /\.fa/ || $output =~ /\.fasta/) {
			print "Err: -o should be FASTQ format (.fq or .fastq) or FASTA format (.fa or .fasta)";
			printHelp();	
		}
	}
	open (OFH, ">$output") or die $!;
	$to_prnt = 1;
}

unless ($ihit && $ifqlist) {
	print "\nErr: -ihit and -ifqlist should be specified\n\n";
	printHelp();
	exit;
}

my $min_len = 1000000000;
my $max_len = 0;
my $num_reads = 0;
my $sum_len;

# ------------- read through $ihit file -------------
unless ($option{silent}) {
	print "\tAnalyze hit file: $ihit\n";
}
my %rname_hash;
open (HIT, $ihit) or die $!;
while (my $line = <HIT>) {
	next if ($line =~ /^\s*$/); #skip blank lines
	chomp($line);
	my @entries = split ('\t', $line);
	$rname_hash{$entries[0]} = 1;	
	++ $num_reads;
}
close (HIT);
print "\n\tnum_reads in $ihit: $num_reads\n\n";
$num_reads = 0;

# ---------- search in fastq files ------------------
my %qualcnts;
foreach my $ifq (@fqfiles) {
	unless ($option{silent}) {
		print "\tParse fastq file: $ifq\n";
	}
	
	open (FQ, $ifq) or die $!;
	my $idx = 0;
	my $to_analyze = 0;
	while (my $line = <FQ>) {
		chomp ($line);
		if ($idx == 0) {
			if ($line =~ /^\@/) { 
				$line = substr($line, 1, length($line) - 1); # remove the first char of readname to be consistent with bowtie output 
				if (exists $rname_hash{$line}) {
					$to_analyze = 1;
					++ $num_reads;
					if ($to_prnt) {
						if ($is_fastq) { print OFH "\@$line\n"; }
						else { print OFH ">$line\n"; }
					}
				} 
				$idx = 1;
			} else { 
				die "\npotential error in the format of file $ifq\n\n";	
			}	
		} elsif ($idx == 1) {
			if ($to_analyze == 1) {
				my $len = length ($line);
				if ($len < $min_len) { $min_len = $len; }
				if ($len > $max_len) { $max_len = $len; }
				$sum_len += $len;
				if ($to_prnt) {
					if ($is_fastq) { print OFH "$line\n"; }
					else { print OFH "$line\n"; }
				}				
			}
			$idx = 2;
		} elsif ($idx == 2) {
			$idx = 3;
			if ($to_analyze) {
				if ($to_prnt) {
						if ($is_fastq) { print OFH "\@$line\n"; }
				}
			}
		} elsif ($idx == 3) { # quality scores
			if ($to_analyze == 1) {
				my $len = length ($line);
				for (my $i = 0; $i < $len; ++ $i) {
					my $base = substr ($line, $i, 1);
					$base = ord($base) - 33;
					if (exists $qualcnts{$base}) { ++ $qualcnts{$base}; }
					else { $qualcnts{$base} = 1; }
				}
				if ($to_prnt) {
					if ($is_fastq) { print OFH "\@$line\n"; }
				}			
			}
			$to_analyze = 0;
			$idx = 0;
		}
	}
	close (FQ);
}

if ($output) { close (OFH); }

my $avg = sprintf("%0.3f", ($sum_len/$num_reads));
if ($num_reads == 0) { print "No reads found\n"; exit;}
print "\n\tNum Reads : $num_reads\n";
print "\tAverage : ".$avg."\n";
print "\tMinimum Length : ".$min_len."\n";
print "\tMaximum Length : ".$max_len."\n";
print "\n\tQuality Score Histogram\n";
foreach my $key  (sort {$b <=> $a} keys %qualcnts) {
	print "\t".$key."\t".$qualcnts{$key}."\n";
}

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief:	Given read names specified in a tsv file as the 1st field and a list of fastq files\n\t
		  	separated by comma in the commandline, searching the reads (optional output as FASTA\n\t
		  	or FASTQ format if -o) in the fastq files and	generate the read statistics.\n\t
		  	Note: the query read name should not contain the first @ char of fq format.\n\n";
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

 