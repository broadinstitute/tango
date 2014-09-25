# This program automatically downloads specified genomes from NCBI genome ftp database 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h     		    => '',
	quiet			=> '', 
	refdir			=> '',  # I/O setting
	taxdir			=> '', 
);

GetOptions(
	"h"					=> \$option{h},
	"quiet"				=> \$option{quiet},
	"refdir=s"			=> \$option{refdir},
	"taxdir=s"			=> \$option{taxdir},
) || printHelp (); 

if ($option{h}) { printHelp();}
my $refdir = $option{refdir};
my $taxdir = $option{taxdir};

unless ($refdir && $taxdir) {
	print "\nErr: -refdir and -taxdir should be specified\n";
	printHelp ();
}


my $cmd_get_bacteria = "perl /seq/viral/analysis/xyang/FUO/scripts/fetch_ncbiftp_bacteria.pl";
my $cmd_get_human = "perl /seq/viral/analysis/xyang/FUO/scripts/fetch_ncbiftp_human.pl";
my $cmd_get_virus = "perl /seq/viral/analysis/xyang/FUO/scripts/fetch_ncbiftp_viruses.pl";
my $cmd_get_fungi = "perl /seq/viral/analysis/xyang/FUO/scripts/fetch_ncbiftp_fungi.pl";
my $cmd_get_taxonomy = "perl /seq/viral/analysis/xyang/FUO/scripts/fetch_ncbi_taxonomy.pl";
my $cmd_clean_taxa_name = "perl /seq/viral/analysis/xyang/FUO/scripts/clean_taxa_names.pl";
my $cmd_prepdb = "/seq/viral/analysis/xyang/FUO/scripts/PrepDB/bin/prepdb"; 
#my $cmd_skeletondb = "/seq/viral/analysis/xyang/FUO/scripts/SkeletonDB/bin/skeletondb";

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = $mon."-".$mday."-".$year;

if (-d $taxdir) { `rm -i -r $taxdir`; }
`mkdir $taxdir`;

system ("$cmd_get_taxonomy -odir $taxdir");
system ("$cmd_clean_taxa_name -i $taxdir/names.dmp -o $taxdir/names.clean.dmp");

unless ($option{quiet}) { print_time ("Download NCBI refs"); }

if (-d $refdir) { `rm -i -r $refdir`; }
`mkdir $refdir`;

system ("$cmd_get_human -odir $refdir");
system ("$cmd_get_virus -odir $refdir");
system ("$cmd_get_fungi -odir $refdir");
system ("$cmd_get_bacteria -odir $refdir");
system ("$cmd_get_bacteria -odir $refdir -draft");

# move protein & mrna to different folder 
opendir (DIR, $refdir) or die $!;
my $rawdir = "$refdir/raw";
my $prdir = "$refdir/raw/pr";
my $mrnadir = "$refdir/raw/mrna";
my $genomedir = "$refdir/raw/genome";
my $plasdir = "$refdir/raw/plasmid";

# generate dir for protein & mRNA
if (-d "$rawdir") { `rm -rf $rawdir`; }
`mkdir $rawdir` == 0 or die $!;
`mkdir $prdir` == 0 or die $!;
`mkdir $mrnadir` == 0 or die $!;
`mkdir $genomedir` == 0 or die $!;
`mkdir $plasdir` == 0 or die $!;

while (my $file = readdir(DIR)) {
	if ($file =~ /\.pr$/) { `mv $refdir/$file $prdir`; }
	elsif ($file =~ /-rna/) { `mv $refdir/$file $mrnadir`; }
	elsif ($file =~ /-plasmid/) { `mv $refdir/$file $plasdir`; }
	elsif ($file =~ /-genome/) { `mv $refdir/$file $genomedir`; } 
}

# invoke prepdb and spilt db by species level 
unless ($option{quiet}) { print_time ("PrepDB"); }

print "[CMD] $cmd_prepdb -iginode $taxdir/gi_taxid_nucl.dmp -oginode $taxdir/gi_range_taxid.txt\n";
system ("$cmd_prepdb -iginode $taxdir/gi_taxid_nucl.dmp -oginode $taxdir/gi_range_taxid.txt");


my $species_genome_dir = "$refdir/species-genome/";
if (-d $species_genome_dir) { `rm -rf $species_genome_dir`; } 
 
`mkdir $species_genome_dir`;
unless ($option{quiet}) { print_time ("Split genome sequences at the Species level"); }
print "[CMD] $cmd_prepdb -iginoderaw $taxdir/gi_taxid_nucl.dmp -iginode $taxdir/gi_range_taxid.txt -itree $taxdir/nodes.dmp -idbdir $genomedir -odbdir $species_genome_dir\n";
system ("$cmd_prepdb -iginoderaw $taxdir/gi_taxid_nucl.dmp -iginode $taxdir/gi_range_taxid.txt -itree $taxdir/nodes.dmp -idbdir $genomedir -odbdir $species_genome_dir");

my $species_mrna_dir = "$refdir/species-mrna/";
if (-d $species_mrna_dir) {	`rm -rf $species_mrna_dir`; }

`mkdir $species_mrna_dir`;
unless ($option{quiet}) { print_time ("Split mRNA sequences at the Species level"); }
print "[CMD] $cmd_prepdb -iginoderaw $taxdir/gi_taxid_nucl.dmp -iginode $taxdir/gi_range_taxid.txt -itree $taxdir/nodes.dmp -idbdir $mrnadir -odbdir $species_mrna_dir\n";
system ("$cmd_prepdb -iginoderaw $taxdir/gi_taxid_nucl.dmp -iginode $taxdir/gi_range_taxid.txt -itree $taxdir/nodes.dmp -idbdir $mrnadir -odbdir $species_mrna_dir");


unless ($option{quiet}) { print_time ("DB update complete"); }


################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Update database: retrieve references and taxonomy from NCBI\n";
		print "\nusage: ./$0 -refdir [refdir] -taxdir [taxdir]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-refdir: output DIR for references, NOTE this dir will firstly be cleaned\n";
		print "-taxdir: output DIR for taxonomy, NOTE this dir will firstly be cleaned\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

sub print_time {
	my $msg = shift;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	print "\n$msg - $hour:$min\n\n"; 
}