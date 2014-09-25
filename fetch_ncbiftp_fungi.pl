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
	"h"							=> \$option{h},
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
unless ($option{quiet}) { print "\nDownload Fungi - $hour:$min\n"; }

# ---- generate temporary director ----------
my $otmpDIR = $odir."/ncbi-ftp-fungi";
if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";

my $ftplink = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/";

# ---- get genomes, plasmid, rrna, mt ----------------
my $type = "genomic.fna.gz"; 
my $output = $odir."/fungi-genome-$print_date.fa";

combineftpgz ($ftplink, $type, $output, $otmpDIR);

unless ($option{quiet}) {
	system("grep \">\" $output | wc -l");
	print "\tgenomic sequences recorded to $output\n";	
}

# ---- protein ---- 
$output = $odir."/fungi-$print_date.pr"; 
$type = "protein.faa.gz";  #protein
combineftpgz ($ftplink, $type, $output, $otmpDIR);

unless ($option{quiet}) {
	system("grep \">\" $output | wc -l");
	print "\tprotein sequences recorded to $output\n";	
}

# ---- RNA ----
$output = $odir."/fungi-rna-$print_date.fa"; 
$type = "rna.fna.gz"; # mRNA
combineftpgz ($ftplink, $type, $output, $otmpDIR);
unless ($option{quiet}) {
	system("grep \">\" $output | wc -l");
	print "\tmRNA sequences recorded to $output\n";	
}
unless ($option{quiet}) { print "\nDONE! \n";}

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve Fungal genomes, plasmids, mRNAs, and Mt from NCBI database\n";
		print "\nusage: ./$0 -odir [odir]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbiftp_fungi.pl -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

sub findFiles {
  my $filetype = ${$_[0]}{ft};
  my $ref_files = ${$_[1]}{fs};
  
	my $file = $File::Find::name;	
	if ($file =~ /\.($filetype)$/i) {
#			push @files, $file;
		push @ {$ref_files}, $file;
	}
}

# obtain gz files from a ftp folder and unzip to output file
sub combineftpgz {
	my $ftplink = shift;
	my $type = shift;
	my $output = shift;
	my $otmpDIR = shift;
	
	if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
	if (-e $output) { system ("rm $output"); }

	system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";
	
	system ("wget -q -A \".$type\" -m -nd $ftplink -P $otmpDIR");
	
	my @target_files;
	find( sub { findFiles({ ft => $type }, { fs => \@target_files })}, $otmpDIR);
	foreach my $file (@target_files) { system ("gunzip -c $file >> $output"); }
	
	system ("rm -rf $otmpDIR");
}
 