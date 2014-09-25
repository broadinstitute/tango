# @description: given a tsv file specifying read names in the first column, 
# remove from an input fastq file any read matching input names.
# The remaining reads are written to an output file whereas the removed
# fq reads are ADDED to an optional output.
# Assume in fastq file, each read name starts by @ and ends by /1 or /2, 
# whereas in tsv file, read name contains neither

#!/usr/bin/perl
use strict;
use Getopt::Long;

my %option = (
	h         => '',
	quiet 		=> '',
	# I/O setting 
	itsv			=> '', 
	ifq				=> '',
	ofq				=> '',
	ofqrm 		=> '',
);

GetOptions(
	"h"					=> \$option{h},
	"quiet"				=> \$option{quiet},
	"itsv=s"				=> \$option{itsv},
	"ifq=s"		  		=> \$option{ifq},
	"ofq=s" 			=> \$option{ofq},
	"ofqrm=s"			=> \$option{ofqrm},
) || printHelp (); 

if ($option{h}) { printHelp();}

my $itsv = $option{itsv};
my $ifq = $option{ifq};
my $ofq = $option{ofq};
my $ofqrm = $option{ofqrm};

unless ($itsv && $ifq && $ofq) {
	print "\nErr: -itsv, -ifq and -ofq should be specified\n";
	printHelp ();
}


unless ($option{quiet}) { 
	print "\n[CMD]\tperl rm_read_by_name.pl -itsv $itsv -ifq $ifq -ofq $ofq\n\n"; 
}

my $num_reads = 0;

# ------------- read through $itsv file -------------
unless ($option{quiet}) {	print "\n\tAnalyze itsv file: $itsv\n\n"; }
my %rname_hash;
open (ITSV, $itsv) or die $!;
while (my $line = <ITSV>) {
	next if ($line =~ /^\s*$/); #skip blank lines
	chomp($line);
	my @entries = split ('\t', $line);
	my $rname = $entries[0];
	# rm first char '@' and end /1 or /2 if applicable 
	if (substr($rname, 0, 1) eq '@') { $rname = substr($rname, 1, length($rname) - 1); }
	if ($rname =~ /\/1$/ || $rname =~ /\/2$/) {	$rname = substr ($rname, 0, length ($rname) - 2); 	}
	$rname_hash{$rname} = 1;	
	++ $num_reads;
}
close (ITSV);
print "\tnum_reads in $itsv: $num_reads\n\n";

open my $fh_ofq, ">$ofq" or die $!; # output fq file
my $fh_ofqrm;
if ($ofqrm) {
	open $fh_ofqrm, ">>$ofqrm" or die $!; # output fqrm file
}

unless ($option{quiet}) {	print "\tParse input fastq file: $ifq\n\n";	}
open (FQ, $ifq) or die "Cannot open $ifq\n";
while (my $line = <FQ>) {
		my @fq_entry;
		if ($line =~ /^\@/) {
			push (@fq_entry, $line);
			for (my $i = 0; $i < 3; ++ $i) {
				$line = <FQ>;
				push (@fq_entry, $line);
			}
			my $header = $fq_entry[0];
			chomp ($header);
			# rm first char '@'
			if (substr($header, 0, 1) eq '@') { $header = substr($header, 1, length($header) - 1); }
			if ($header =~ /\/1$/ || $header =~ /\/2$/) {
				# rm /1 or /2
				$header = substr ($header, 0, length ($header) - 2); 
			}
			
			if (!exists $rname_hash{$header}) { # not found
				output_read (\@fq_entry, $fh_ofq);
			}	else {
				if ($ofqrm) {
					output_read (\@fq_entry, $fh_ofqrm);
				}
			}
		} else { die "\npotential error in the format of file $ifq\n\n"; }
} # while 
close (FQ);
close ($fh_ofq);
if ($ofqrm) {	close ($fh_ofqrm); }

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: given a tsv file specifying read names in the first column, \n";
		print " remove from an input fastq file any read matching input names.\n";
		print " The remaining reads are written to an output file whereas the removed\n";
		print " fq reads are ADDED to an optional output.\n";
		print " Assume in fastq file, each read name starts by @ and ends by /1 or /2,\n";
		print " whereas in tsv file, read name contains neither\n\n";
		print "\nusage: ./$0 -itsv [names.tsv] -ifq [input.fq] -ofq [output.fq]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-itsv: input file in tsv format with 1st field being read names\n";
		print "-ifq: input fq file\n";
		print "-ofq: output fq file storing remaining reads\n";
		print "-ofqrm: (optional) output fq file storing removed reads\n";		
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

####################################
# output fastq read entry to a fastq file handle
sub output_read {
	my $array_ref = shift;
	my $fh = shift;
	foreach my $elem (@$array_ref) {
		print $fh $elem;
	}
} # output_read