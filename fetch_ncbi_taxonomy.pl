# This program automatically downloads NCBI taxonomy 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h       => '',
	quiet	=> '', 
	odir	=> '',  # I/O setting
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
unless ($option{quiet}) { print "\nDownload NCBI Taxonomy - $hour:$min\n"; }

my $ftplink = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/";

# taxdmp 
my $downloadfile = $ftplink."taxdmp.zip";
system ("wget -q $downloadfile -P $odir");
system ("unzip -q -o $odir/taxdmp.zip -d $odir"); # -o replacing existing file if available
`rm $odir/taxdmp.zip`;

# gi_taxid_nucl

$downloadfile = $ftplink."gi_taxid_nucl.dmp.gz";
system ("wget -q $downloadfile -P $odir");

system ("gunzip -q -d -f $odir/gi_taxid_nucl.dmp.gz"); #-d decompress -f force overwrite

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve Taxonomy files \n";
		print "\nusage: ./$0 -odir [odir]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbiftp_taxonomy.pl -odir /seq/viral/analysis/xyang/FUO/curated_database/taxonomy\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}
