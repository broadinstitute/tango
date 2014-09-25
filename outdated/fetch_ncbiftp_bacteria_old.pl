# This program automatically downloads specified genomes from NCBI genome ftp database 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h         => '',
	quiet			=> '', 
	draft			=> '',
	pr				=> '',
	odir			=> '',  # I/O setting
);

GetOptions(
	"h"						=> \$option{h},
	"quiet"					=> \$option{quiet},
	"draft"					=> \$option{draft},
	"pr"					=> \$option{pr},
	"odir=s"				=> \$option{odir}
) || printHelp (); 

if ($option{h}) { printHelp();}
my $odir = $option{odir};
my $type = "fna";
if ($option{pr}) { $type = "faa"; }

unless ($odir) {
	print "\nErr: -odir should be specified\n";
	printHelp ();
}

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = $mon."-".$mday."-".$year;
unless ($option{quiet}) {	print "\nDownload Bacteria - $hour:$min\n"; }

# ---- generate temporary director ----------
my $otmpDIR = $odir."/ncbi-ftp-genome-Bacteria";

if ($option{draft}) { $otmpDIR .= "-draft"; }
if (-d $otmpDIR) { system ("rm -rf $otmpDIR"); }
system ("mkdir $otmpDIR") == 0 or die "creating $otmpDIR failed\n";

my $ftplink = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria";

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
	
	unless ($option{quiet}) {	print "\n".@tgz_files." tgz files to be processed\n";	}
	
	my $cnt = 0; 
	foreach my $file (@tgz_files) {
		system ("tar -xzf $file -C $otmpDIR");
		system ("rm $file");
		++ $cnt;
		if ($cnt % 200 == 0) { unless ($option{quiet}) { print "\t$cnt files processed\n"; } } 
	}
 
	unless ($option{quiet}) {
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
		print "\nTime - $hour:$min\n";
	}

} else {
	# ---- get genomes and unzip ----------------
	
	my $ofile = $otmpDIR."/all.$type.tar.gz";
	my $link = $ftplink."/all.$type.tar.gz";	
	get_extract_link ($link, $ofile, $otmpDIR);
	unless ($option{quiet}) {
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
		print "\nTime - $hour:$min\n";
	}
} # else

my $oAll = $odir."/bacteria-all-$print_date.fa";
if ($option{draft}) { $oAll = $odir."/bacteria-all-draft-$print_date.fa"; }

if (-e $oAll) { `rm $oAll`; }

print "\nfind $otmpDIR -name \"*.$type\" -exec cat {} \+ > $oAll\n";
system ("find $otmpDIR -name \"*.$type\" -exec cat {} \+ > $oAll");
#my @files;
#find( sub { findFiles({ ft => $type }, { fs => \@files })}, $otmpDIR);
#unless ($option{quiet}) {	print "\nmerge ".@files." files to $output\n"; }
#foreach my $file (@files) { 
#	system ("cat $file >> $output"); 
#	`rm $file`;
#}

# ---- split file into Plasmids, genome, RNA -----
print "\nSplit $oAll into genomic and plasmids\n";
my $num_genomes = 0;
my $num_plasmids = 0;

my $oGenome = $odir."/bacteria-genome-$print_date.fa";
my $oPlasmid = $odir."/bacteria-plasmid-$print_date.fa";
if ($option{draft}) { 
	$oGenome = $odir."/bacteria-genome-draft-$print_date.fa";
	$oPlasmid = $odir."/bacteria-plasmid-draft-$print_date.fa";
}

open(my $ofh_gen, '>', $oGenome) or die "Could not open file $oGenome to write\n";
open(my $ofh_plas, '>', $oPlasmid) or die "Could not open file $oPlasmid to write\n";

open (my $fh_all, $oAll) or die "cannot read $oAll\n";	

my $header;
my $seq = "";
while (my $line = <$fh_all>) {
	next if ($line =~ /^\s*$/); #skip blank lines

	if ($line =~ /\>/) {
		if ($header =~ /plasmid/i) {
			++ $num_plasmids;
			if ($seq ne "") { 
				print $ofh_plas $header;
				print $ofh_plas $seq;	
			}
		}	else { 
			if ($seq ne "") { 
				print $ofh_gen $header;
				print $ofh_gen $seq;
				++ $num_genomes;
			}			
		} 
		$header = $line;
		$seq = "";
	} else { $seq .= $line; }
}

# the last sequence
if ($header =~ /plasmid/i) {
	++ $num_plasmids;
	if ($seq ne "") { print $ofh_plas $header.$seq;	}
}	else { 
	if ($seq ne "") { 
		print $ofh_gen $header.$seq; 
		++ $num_genomes;
	}			
} 

my $num_total = $num_genomes + $num_plasmids;
print "\nNum of Genomic, Plasmids: $num_genomes, $num_plasmids ($\n $num_total)\n";		
close ($ofh_gen);
close ($ofh_plas);
close ($fh_all);


unless ($option{quiet}) {	print "\nclean up: rm -rf $otmpDIR\nrm $oAll";}
system ("rm -rf $otmpDIR");
system ("rm $oAll");
unless ($option{quiet}) { print "\nDONE! \n";}
exit;

### the following part is disabled right now, which go through each file and split
# mt, genome, plasmid, and rna into different files.

# ---- Create output genome /plasmid /RNA /Mitochondria files ----
my $oGenome = $odir."/bacteria-genome-$print_date";
my $oPlasmid = $odir."/bacteria-plasmid-$print_date";
my $oMitochondria = $odir."/bacteria-mitochondria-$print_date";
my $oRNA = $odir."/bacteria-rna-$print_date";
if ($option{draft}) {
	$oGenome .= "-draft";
	$oPlasmid .= "-draft";
	$oMitochondria .= "-draft";
	$oRNA .= "-draft";
}

if ($option{pr}) {
	$oGenome .= "-pr.fa";
	$oPlasmid .= "-pr.fa";
	$oMitochondria .= "-pr.fa";
	$oRNA .= "-pr.fa";
} else {
	$oGenome .= ".fa";
	$oPlasmid .= ".fa";
	$oMitochondria .= ".fa";
	$oRNA .= ".fa";
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
				print "\n\tTime - $hour:$min\n";
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
	unless ($option{quiet}) {	print "\n$num_genomes genomes written to $oGenome\n";	}
}
if ($num_plasmids == 0) {	system ("rm $oPlasmid"); }
else {	unless ($option{quiet}) {	print "\n$num_plasmids plasmids written to $oPlasmid\n"; }
}
if ($num_rna == 0) { system ("rm $oRNA"); }
else {
	unless ($option{quiet}) {	print "\n$num_rna RNAs written to $oRNA\n";	}
}
if ($num_mito == 0) { system ("rm $oMitochondria"); }
else {
	unless ($option{quiet}) {	print "\n$num_mito mitochondria written to $oMitochondria\n";	}
}

unless ($option{quiet}) {	print "\nclean up: rm -rf $otmpDIR\n\n";}
system ("rm -rf $otmpDIR");

unless ($option{quiet}) { print "\nDONE! \n";}

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve bacteria genomes and plasmids from NCBI database\n";
		print "\nusage: ./$0 -odir [odir]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-draft: if specified, draft genomes will be retrieved\n";
		print "-pr: by default nt is retrieved but if specified, pr will be retrieved\n";
		print "-odir: output DIR\n";
		print "\nExample: perl fetch_ncbiftp_bacteria.pl -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

sub findFiles {
  my $filetype = ${$_[0]}{ft};
  my $ref_files = ${$_[1]}{fs};
  
	my $file = $File::Find::name;	
	if ($file =~ /\.($filetype)$/i) {
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

