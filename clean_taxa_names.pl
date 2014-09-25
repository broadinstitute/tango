#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;
use File::stat;
use File::Find;
	
my %option = (
	h  		       => '',
	quiet			=> '', 
	i				=> '',  # I/O setting
	o				=> '',
);

GetOptions(
	"h"						=> \$option{h},
	"quiet"					=> \$option{quiet},
	"i=s"						=> \$option{i},
	"o=s"					=> \$option{o},
) || printHelp (); 

my $in = $option{i};
my $out = $option{o};

unless ($in && $out) { print "\n[ERR] both -i and -o should be specified\n"; printHelp();	exit; }

unless ($option{quiet}) { 
	print "\n[CMD]\tperl $0 -i $in -o $out\n\n"; 
}

open(INPUT, $in) or die "unable to open $in to read\n\n";
open(OUTPUT, ">$out") or die "can't open $out to write\n";	
	
my $prev_id = -1;
my $is_scientific_name_found = 0;
my @prev_id_entry;

while (<INPUT>) {
	next if /^(\s)*$/;  # skip blank lines
	my $line = $_;
	unless ($line =~ /misspelling/ || $line =~ /type material/ || $line =~ /authority/) {		
		my @fields = split (/\s/, $line);
		my $id = $fields[0];
		if ($prev_id == -1) { # the first entry
			$prev_id = $id; 
			if ($line =~ /scientific name/) {	
				$is_scientific_name_found = 1; 
				print OUTPUT $line;
			} else {	@prev_id_entry = ($line); }			
		} else {
			if ($id == $prev_id) {			
				unless ($is_scientific_name_found == 1) {
					if ($line =~ /scientific name/) {	
						$is_scientific_name_found = 1; 
						print OUTPUT $line;		
						@prev_id_entry = ();				
					} else { push (@prev_id_entry, $line); }
				}
			} else {
				foreach (@prev_id_entry) { print OUTPUT $_; }
				$is_scientific_name_found = 0; 
				$prev_id = $id;
				if ($line =~ /scientific name/) {	
					$is_scientific_name_found = 1; 
					print OUTPUT $line;
				} else { @prev_id_entry = ($line); }			
			}
		} # else
	} # unless 
} # while 

# process the last entry
foreach (@prev_id_entry) { print OUTPUT $_; }

close (INPUT);
close (OUTPUT);

unless ($option{quite}) { print "Done!\n\n"; }

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: keep only scientific name if applicable and remove non-informative entries \n(containing certain words e.g. misspelling) in the taxaonomic name file names.dmp\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-i: input file\n";
		print "-o: output file\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

 
 