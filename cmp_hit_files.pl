# @brief
# description: given a folder storing a set of tsv files (these files should have string 'hit'
# in the file name) in the format of [rname][tab][strand][tab][reference][tab][aln_pos]
# Identify pair-wisely the common read names (optional: output to .shared.txt if -prnt), 
#	summarize hits to each reference (output to .refcnts.txt)

#!/usr/bin/perl

use strict;
use Getopt::Long;


#--------------------------------------------------------------------------------------------------

my %option = (
	h       => '',
	silent 	=> '',
	# I/O setting 
	idir		=> '',
	prnt		=> '', 
);

GetOptions(
	"h"				=> \$option{h},
	"silent" 	=> \$option{silent},
	"idir=s"	=> \$option{idir}, 
	"prnt"		=> \$option{prnt}, 
) || printHelp (); 

if ($option{h}) { printHelp();}

my $idir = $option{idir};
my $silent_option = "";
if ($option{silent}) { $silent_option = " -silent "; }
unless ($idir) {	print "\n\tERR:  -idir should be specified\n"; printHelp ();}

# identify all bam files 
opendir (DIR, $idir) or die $!;	
# first get all files 
my @file_names;
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./); # ignore . .. files
	# obtain prefix of the file
	if ($file =~ /hit/) { push (@file_names, $file); }	
}
close (DIR);

my @array_of_hash;
my $num_files = $#file_names + 1;
for (my $i  = 0; $i < $num_files; ++ $i) {
	my %rname_hash;
	my %ref_rcnt; # refname -> num of read hits
	my %ref_pos_cnt; # refname -> aln_pos -> cnt;
	my $path = $idir."/$file_names[$i]";
	open (FH, $path) or die $!;
	while (my $line = <FH>) {
		chomp($line);
		my @entries = split ('\t', $line); # [readname] [strand] [refname] [aln_pos]
		$rname_hash{$entries[0]} = 1;
		if (exists $ref_rcnt{$entries[2]}) { ++ $ref_rcnt{$entries[2]}; }
		else { $ref_rcnt{$entries[2]} = 1; }
		if (exists $ref_pos_cnt{$entries[2]}{$entries[3]}) {  ++ $ref_pos_cnt{$entries[2]}{$entries[3]}; }
		else {  $ref_pos_cnt{$entries[2]}{$entries[3]} = 1; }
	}
	close (FH);
	push (@array_of_hash, \%rname_hash);
	
	#output %ref_hash to file 
	$file_names[$i] =~ /(.+)\.(.+)/;
	my $output_file = $idir."/$1".".refcnts.txt";
	$output_file =~ s/\.hit//;
	output_ref_hash (\%ref_rcnt, \%ref_pos_cnt, $output_file);
}

print "list of files analyzed: @file_names\n";
closedir DIR;

for (my $i  = 0; $i < $num_files; ++ $i) {
	my $base;
	if ($option{prnt}) {
		$file_names[$i] =~ /(.+)\.(.+)/;
		$base = $idir."/$1";
	}
	for (my $j = $i + 1; $j < $num_files; ++ $j) {
		my $ofile;
		if ($option{prnt}) {
			$file_names[$j] =~ /(.+)\.(.+)/;
			$ofile = $base."."."$1".".shared.txt";
			print "ofile = $ofile\n";
			$ofile =~ s/\.hit//g;
			print "ofile = $ofile\n";
		}	
		my $num_comm_elems = get_comm_elements ($array_of_hash[$i], $array_of_hash[$j], $ofile); 
		print "$file_names[$i]\t$file_names[$j]\t$num_comm_elems\n";
		
	}
}
#################################################################################################################################
# sub-functions 
################################
pairwise comparing the hit files, obtain shared elements (.shared.txt), summarize hits to each reference (.refcnts.txt) 
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: given a folder storing a set of tsv files (these files should have string 'hit'\n";
		print "\tin the file name) in the format of [rname][tab][strand][tab][reference][tab][aln_pos]\n";
		print "\tIdentify pair-wisely the common read names (optional: output to .shared.txt if -prnt)\n";
		print "\tsummarize hits to each reference (output to .refcnts.txt). Note: the output read names\n";
		print "\tdo not contain the first @ char of fq format. \n\n";
		
		print "usage: ./$0 {-idir [inDIR]} \n\n";
		print "-h: help; for cmd options\n";
		print "-silent: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-idir: input DIR containing bowtie result files (file name contain the string 'hit') in the format:\n";
		print "\t [rname][tab][strand][tab][reference][tab][aln_pos]\n";
		print "-prnt: opt; if specified, output common read(names) to file\n"; 
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

sub get_comm_elements {
	my $hash_ref = shift;
	my $hash_ref2 = shift;
	my $ofile = shift;
	if ($ofile) {
		open (FH, ">$ofile") or die $!;
	}
	my $cnt = 0;
	foreach my $key (keys %$hash_ref) {
		if (exists $hash_ref2->{$key}) {
			++ $cnt;
			if ($ofile) {	
				print FH $key."\n";
			}
		}
	}	
	if ($ofile) { close (FH); }
	return $cnt;
}

sub output_ref_hash {
	my $reference_hash_refcnt = shift;
	my $reference_hash_refposcnt = shift;	
	my $ofile = shift;
	open (FH, ">$ofile") or die $!;
	foreach (sort { ($reference_hash_refcnt->{$b} <=> $reference_hash_refcnt->{$a}) } keys %$reference_hash_refcnt) {
		my $ref = $_;
		print FH "$ref\t$reference_hash_refcnt->{$ref}\n";
		foreach my $key (sort keys %{$reference_hash_refposcnt->{$ref}}) {
			print FH "$key\t$reference_hash_refposcnt->{$ref}{$key}\n"; 
		} 
	}	
	close (FH);
}