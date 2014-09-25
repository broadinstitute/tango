# This program automatically downloads Human genome, transcripts, MT from NCBI ftp database 

# 1. Human genome ref including a) RefSeq, b) MT, c) unlocalized seq, d) unplaced seq, e) alternate locus group
# 2. Transcripts ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/README -- "The RNA and protein directories provide sequence files in three
#    formats representing all of the mRNA, non-coding transcript, and protein model reference sequences (RefSeq) exported as part of the
#			genome annotation process.

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h         => '',
	quiet			=> '', 
	chr_pattn => "hs_ref_GRCh*_*.fa.gz",
	odir			=> '',  # I/O setting
);

GetOptions(
	"h"							=> \$option{h},
	"quiet"					=> \$option{quiet},
	"chr_pattn" 		=> \$option{chr_pattn},
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
unless ($option{quiet}) {	print "\nDownload Human - $hour:$min\n"; }

# ---- generate temporary director ----------
my $otmpDIR = $odir."/ncbi-ftp-human";
if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";

# ---- 1. human genome assembly (chromosomes, MT, alternate locus group, unlocalized and unplaced seq --------
my $ftp_link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/";
my $chr_pattn = $option{chr_pattn};
my $output = $odir."/human-genome-$print_date.fa";

#unless($option{quiet}) { print "\n\twget -q -A \"$chr_pattn\" -m -nd $ftp_link -P $otmpDIR\n";	}
system ("wget -q -A \"$chr_pattn\" -m -nd $ftp_link -P $otmpDIR");

my @target_files;
find( sub { findFiles({ ft => "fa.gz" }, { fs => \@target_files })}, $otmpDIR);
if (-e $output) { system ("rm $output"); }
#unless ($option{quiet}) { print "\n\t".@target_files." files identified, combined to $output\n"; }
foreach my $file (@target_files) {	system ("gunzip -c $file >> $output"); }

unless ($option{quiet}) {
	system("grep \">\" $output | wc -l");
	print "\tgenomic sequences downloaded\n";	
}

# ---- 2. Annotated transcripts mRNA + non-coding ---- 
$output = $odir."/human-rna-$print_date.fa";
$ftp_link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/RNA/";

combineftpgz ($ftp_link, "fa.gz", $output, $otmpDIR);

unless ($option{quiet}) {
	system("grep \">\" $output | wc -l");
	print "\tmRNA sequences downloaded\n";	
}

unless ($option{quiet}) { print "\nDONE! \n";}

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
		print "-chr_pattn: by default hs_ref_GRCh*_*.fa.gz files are retrieved\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbiftp_human.pl -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
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