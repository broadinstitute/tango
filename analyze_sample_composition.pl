# This program identify all the taxa_rank_SAMPLENAME.txt files in a given directory
# and analyze composition across these files  

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

use File::Find;
	
my %option = (
	h         => '',
	quiet	  => '', 
	idir	  => '',  # I/O setting
	pattern	  => '', 
	ofile	  => '', 
);

GetOptions(
	"h"			=> \$option{h},
	"quiet"		=> \$option{quiet},
	"idir=s"	=> \$option{idir},
	"pattern=s" => \$option{pattern},
	"ofile=s" 	=> \$option{ofile}
) || printHelp (); 

if ($option{h}) { printHelp();}
my $idir = $option{idir};
my $ofile = $option{ofile};
my $pattern = $option{pattern};
unless ($idir && $ofile && $pattern) {
	print "\nErr: -idir, -pattern and -ofile should be specified\n";
	printHelp ();
}

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = $mon."-".$mday."-".$year;

my @target_files;
find( sub { findFiles({ ft => $pattern }, { fs => \@target_files })}, $idir);

my %taxaID_2_sampleNames;
my %taxaID_2_name;
my @sample_names;
my $sample_index = 0;
foreach my $file (@target_files) {
	my @tmp =  split (/$pattern/, $file);
	my @tmp2 = split (/\./, $tmp[-1]);
	my $sample_name = $tmp2[0];
	push @sample_names, $sample_name;
	
	# analyze the current file 	
	
	open INPUT, $file or die $!;
	while (my $line = <INPUT>) {
		chomp($line);
		#print $line."\n";
		my @entries = split (/\t:\t/, $line);
	
		if ($#entries + 1 == 6) {
			my $taxaID = $entries[1];	
			my $taxaName = $entries[5];
			$taxaID_2_name{$taxaID} = $taxaName;
			
			push(@{$taxaID_2_sampleNames{$taxaID}}, $sample_index);
		}
	}
	close (INPUT);
	$sample_index ++;
}
	
open(my $fh, '>', $ofile) or die "Could not open file $ofile";
foreach my $taxaID (keys %taxaID_2_sampleNames) {
    print $fh $taxaID_2_name{$taxaID}."\t";
    
    foreach (@{$taxaID_2_sampleNames{$taxaID}}) {
    	my $sample_index = $_;
        print $fh $sample_names[$sample_index]."\t";
    }
    print $fh "\n";
}	
close $fh;

unless ($option{quiet}) { print "\nDONE! \n";}

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: identify all the files whose name has a specified pattern (substring) in a given directory\n";
		print "\tand analyze composition across these files\n";
		print "\nusage: perl $0 -idir [idir] -pattern [string] -ofile [ofile]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-idir: input DIR\n";
		print "-pattern: input pattern string\n";
		print "-ofile: output file\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

sub findFiles {
  my $pattern = ${$_[0]}{ft};
  my $ref_files = ${$_[1]}{fs};
  
	my $file = $File::Find::name;	
	if ($file =~ /$pattern/i) {
		push @ {$ref_files}, $file;
	}
}

