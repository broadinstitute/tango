# This program automatically downloads specified genomes from NCBI genome ftp database 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h         => '',
	quiet			=> '', 
	odir			=> '',  # I/O setting
);

GetOptions(
	"h"						=> \$option{h},
	"quiet"					=> \$option{quiet},
	"odir=s"				=> \$option{odir}
) || printHelp (); 

if ($option{h}) { printHelp();}
my $odir = $option{odir};

unless ($odir) {
	print "\nErr: -odir should be specified\n";
	printHelp ();
}

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = $mon."-".$mday."-".$year;
unless ($option{quiet}) {	print "\nDownload Viruses - $hour:$min\n"; }

# ---- generate temporary director ----------
my $otmpDIR = $odir."/ncbi-ftp-refseq-viruses";

if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";

my $downloadfile = "viral.1.1.genomic.fna.gz";
my $output = $odir."/viruses-genome-$print_date.fa";
my $ftplink = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/$downloadfile";

#print "\nwget -q -m -nd $ftplink -P $otmpDIR\n";
system ("wget -q -m -nd $ftplink -P $otmpDIR");
system ("gunzip -c $otmpDIR/$downloadfile > $output");
unless ($option{quiet}) {	
	system("grep \">\" $output | wc -l");
	print " genomic sequences recorded to $output\n"; 
}

$downloadfile = "viral.1.protein.faa.gz"; 	
$output = $odir."/viruses-$print_date.pr";
my $ftplink = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/$downloadfile";
system ("wget -q -m -nd $ftplink -P $otmpDIR");
system ("gunzip -c $otmpDIR/$downloadfile > $output");

unless ($option{quiet}) {	
	system("grep \">\" $output | wc -l");
	print " protein sequences recorded to $output\n";
}

system ("rm -rf $otmpDIR");

unless ($option{quiet}) { print "\nDONE!\n";}


################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve genomes and plasmids from NCBI database\n";
		print "\nusage: ./$0 -odir [odir]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-pr: by default nt is retrieved but if specified, pr will be retrieved\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbiftp_viruses.pl -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}


