# This program automatically downloads specified genomes from NCBI genome ftp database 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h         => '',
	quiet			=> '', 
	orgn			=> '',
	draft			=> '',
	pr				=> '',
	odir			=> '',  # I/O setting
);

GetOptions(
	"h"							=> \$option{h},
	"quiet"					=> \$option{quiet},
	"orgn=s"				=> \$option{orgn},
	"draft"					=> \$option{draft},
	"pr"						=> \$option{pr},
	"odir=s"				=> \$option{odir}
) || printHelp (); 

if ($option{h}) { printHelp();}
my $odir = $option{odir};
my $orgn = $option{orgn};
my $type = "fna";
my $type_ffn = "";  # additional genomes for virus
if ($option{pr}) { $type = "faa"; }

unless ($odir) {
	print "\nErr: -odir should be specified\n";
	printHelp ();
}

unless (($orgn eq "Bacteria") || ($orgn eq "Fungi") || ($orgn eq "Viruses")) {
	print "\nErr: -orgn should be specified as from {Bacteria, Fungi, Viruses}\n";
	exit;
}

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = $mon."-".$mday."-".$year;
unless ($option{quiet}) {	print "\n\tTime: $hour hr\t $min min\n"; }

# ---- generate temporary director ----------
my $otmpDIR = $odir."/genome-".$orgn;

# ---- the following commented are for viral download from genome db, which is currently incomplete in NCBI---
#my $otmpDIR_virus = "";
#if ($orgn eq "Viruses") {
#	$type_ffn = "ffn";  # additional genomes may have name ffn if fna was not found
#	$otmpDIR_virus = $otmpDIR."-ffn";
#}

if ($option{draft}) { $otmpDIR .= "-draft"; }
if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";

#if ($otmpDIR_virus) {	
#	if (-d $otmpDIR_virus) { system ("rm -rf $otmpDIR_virus"); }
#	system ("mkdir $otmpDIR_virus") == 0 or die "creating $otmpDIR_virus failed\n";	
#}

my $ftplink = "ftp://ftp.ncbi.nlm.nih.gov/genomes/".$orgn;

if ($orgn eq "Viruses") {
	# download viruses and return ...
	
	$ftplink = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz";
	print "\nwget -q -m -nd $ftplink -P $otmpDIR\n";
  system ("wget -q -m -nd $ftplink -P $otmpDIR");
 	my $oGenome = $odir."/$orgn-genome-$print_date";
	system ("gunzip -c $otmpDIR/viral.1.1.genomic.fna.gz > $oGenome");
	system ("rm -rf $otmpDIR");
	unless ($option{quiet}) {
		print "\nNum genomes: ";	
		system("grep \">\" $oGenome | wc -l");
		print "\n";	
	}
	unless ($option{quiet}) { print "\nDONE!\n\n";}
	exit;
}


if ($option{draft}) { 
	$ftplink .= "_DRAFT"; 
	unless ($option{quiet}) { print "\nretrieve $ftplink/ to $otmpDIR\n"; }
	print "\nwget -q -A \".$type.tgz\" -m -nd $ftplink -P $otmpDIR\n";
	
  system ("wget -q -A \".$type.tgz\" -m -nd $ftplink -P $otmpDIR");
	
	unless ($option{quiet}) {
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
				print "\n\tTime: $hour hr\t $min min\n";
	}

	# unzip all tgz files in $otmpDIR
	print "\nunzip all .tgz files\n";
	my @tgz_files;
	find (sub {findFiles ({ft => "tgz"}, {fs => \@tgz_files})}, $otmpDIR);
	
	unless ($option{quiet}) {
		print "\n".@tgz_files." tgz files to be processed\n";
	}
	
#	exit;
	
	my $cnt = 0; 
	foreach my $file (@tgz_files) {
		system ("tar -xzf $file -C $otmpDIR");
		system ("rm $file");
		++ $cnt;
		if ($cnt % 100 == 0) {
			unless ($option{quiet}) { print "\t$cnt files processed\n"; }
		}
	}
	
	unless ($option{quiet}) {
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
				print "\n\tTime: $hour hr\t $min min\n";
	}

} else {
	# ---- get genomes and unzip ----------------
	
	my $ofile = $otmpDIR."/all.$type.tar.gz";
	my $link = $ftplink."/all.$type.tar.gz";	
	get_extract_link ($link, $ofile, $otmpDIR);
	
	#if ($otmpDIR_virus) {
	#	$ofile = $otmpDIR_virus."/all.$type_ffn.tar.gz";
	#	$link = $ftplink."/all.$type_ffn.tar.gz";
	#	get_extract_link ($link, $ofile, $otmpDIR_virus);
	#}
	
} # else


# ---- Create output genome /plasmid /RNA /Mitochondria files ----

my $oGenome = $odir."/$orgn-genome-$print_date";
my $oPlasmid = $odir."/$orgn-plasmid-$print_date";
my $oMitochondria = $odir."/$orgn-mitochondria-$print_date";
my $oRNA = $odir."/$orgn-rna-$print_date";
if ($option{draft}) {
	$oGenome .= "-draft";
	$oPlasmid .= "-draft";
	$oMitochondria .= "-draft";
	$oRNA .= "-draft";
}
if ($option{pr}) {
	$oGenome .= "-pr.fasta";
	$oPlasmid .= "-pr.fasta";
	$oMitochondria .= "-pr.fasta";
	$oRNA .= "-pr.fasta";
} else {
	$oGenome .= ".fasta";
	$oPlasmid .= ".fasta";
	$oMitochondria .= ".fasta";
	$oRNA .= ".fasta";
}

if (-e $oGenome) { system "rm $oGenome";}
system ("touch $oGenome");
if (-e $oPlasmid) { system "rm $oPlasmid";}
system ("touch $oPlasmid");
if (-e $oRNA) { system "rm $oRNA";}
system ("touch $oRNA");
if (-e $oMitochondria) { system "rm $oMitochondria";}
system ("touch $oMitochondria");

# --- now go through every Fasta file in the subdirectories of $otmpDIR ---
my @files;
find( sub { findFiles({ ft => $type }, { fs => \@files })}, $otmpDIR);

#print "before size = ".@files."\n";

# --- The following if statement is for viruses to be downloaded from ncbi genome db, 
#which is currently not ready. Identify all .ffn entries that do not have .fna and adds to @files --- 
#if ($otmpDIR_virus) {
	# --- get fna file prefix ----
#	my %fna_prefix;
#	my $index = 0;
#	my $target_idx;
#	foreach my $file (@files) {
#		$file =~ /.*?([^\/]+)\.fna$/;
#		my $prefix = $1;		
#		# known wrong placement in ncbi genome viruses 
#		if ($prefix eq "NC_017312") {	$target_idx = $index;	} 
#		else {	++ $index; }
#		$fna_prefix{$prefix} = 1;
#	}
#	splice (@files, $target_idx, 1); # remove this file
	
#	my @ffn_files;
#	find( sub { findFiles({ ft => $type_ffn }, { fs => \@ffn_files })}, $otmpDIR_virus);
			
#	foreach my $file (@ffn_files) {
#		$file =~ /.*?([^\/]+)\.ffn$/;
#		my $prefix = $1;
#		unless (exists $fna_prefix{$prefix}) {
#			unless ($prefix eq "NC_017312") {
#				push @files, $file;
#			} 
#		} #unless 
#	} # foreach 
#} # if ($otmpDIR_virus)

#print "after size = ".@files."\n";
#exit;

print "\nNum of total files to process: ".@files."\n";
my $num_genomes = 0;
my $num_plasmids = 0;
my $num_rna = 0;
my $num_mito = 0;
my $cnt = 0;
foreach my $file (@files) {
		++ $cnt;
		if ($cnt % 5000 == 0) {
			unless ($option{quiet}) { 
				print "\t$cnt files processed\n"; 
				($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
				print "\n\tTime: $hour hr\t $min min\n";
			}
		}
		open INPUT, "$file" or die "cannot open $file\n";	
		while (<INPUT>) {
			my $line = $_;			
			if ($line =~ /\>/) {
				if ($line =~ /plasmid/i) { 
					system ("cat $file >> $oPlasmid"); 
					system ("rm $file");
					++ $num_plasmids;
				}	elsif ($line =~ /\s[tsrm]rna\s/i) {
					system ("cat $file >> $oRNA");
					system ("rm $file");
					++ $num_rna;	
				}	elsif ($line =~/mitochondria\s/i){
					system ("cat $file >> $oMitochondria");
					system ("rm $file");					
					++ $num_mito;						
				}	else { 
					system ("cat $file >> $oGenome"); 
					system ("rm $file");					
					++ $num_genomes;
				} 
				last;
			}
		}
		close (INPUT);
}

if ($num_genomes == 0) { system ("rm $oGenome"); }
else {
	unless ($option{quiet}) {
		print "\n$num_genomes genomes written to $oGenome\n";		
	}
}
if ($num_plasmids == 0) {	system ("rm $oPlasmid"); }
else {
	unless ($option{quiet}) {
		print "\n$num_plasmids plasmids written to $oPlasmid\n";		
	}
}
if ($num_rna == 0) { system ("rm $oRNA"); }
else {
	unless ($option{quiet}) {
		print "\n$num_rna RNAs written to $oRNA\n";		
	}
}
if ($num_mito == 0) { system ("rm $oMitochondria"); }
else {
	unless ($option{quiet}) {
		print "\n$num_mito mitochondria written to $oMitochondria\n";		
	}
}

unless ($option{quiet}) {	print "\nclean up: rm -rf $otmpDIR\n\n";}
system ("rm -rf $otmpDIR");
#if ($otmpDIR_virus) { 
#	unless ($option{quiet}) {	print "\nclean up: rm -rf $otmpDIR_virus\n\n";}
#	system ("rm -rf $otmpDIR_virus"); 
#}

unless ($option{quiet}) { print "\nDONE! \n";}

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve genomes and plasmids from NCBI database\n";
		print "\nusage: ./$0 -odir [odir] -orgn [organism] -draft\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-orgn: organism choicie in {Bacteria,Viruses,Fungi}\n";
		print "-draft: if specified, draft genomes will be retrieved\n";
		print "-pr: by default nt is retrieved but if specified, pr will be retrieved\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbi_genome_ftp.pl -orgn Bacterial -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
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

# brief wget link and extract to the target folder
sub	get_extract_link {
	my $ftplink = shift;
	my $ofile = shift;
	my $odir = shift;
	
	if (-e $ofile) { system ("rm $ofile"); }
	
	unless ($option{quiet}) { print "\nretrieve $ftplink\n"; }
	system ("wget -q $ftplink -O $ofile") == 0 or die "fetch $ftplink failed\n";

	unless ($option{quiet}) { 
		print "\nunzip retrieved link to $odir\n"; 
		print "\ntar -xf $ofile -C $odir\n";
	}
	system ("tar -xf $ofile -C $odir") == 0 or die "tar ... failed\n";
} # get_extract_link

