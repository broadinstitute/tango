#!/usr/bin/perl

use strict;
use Getopt::Long;
use 5.010; 
use Time::HiRes;

#--------------------------------------------------------------------------------------------------
# set dependent script path
my $sam2fastq_jar = "/seq/software/picard/current/bin/SamToFastq.jar";
my $duprm_c = "/seq/viral/analysis/xyang/programs/DupRm/bin/duprm"; 
#--------------------------------------------------------------------------------------------------

my %option = (
	h           => '',
	silent			=> '',
	p 					=> 8,
	ibam				=> '',
	ifq 				=> '',	 # input fastq pair 1
	ifq2				=> '',   # input fastq pair 2
	batch 			=> 500000, # batch of read pairs to be processed 
	odir				=> '',
	oprefix			=> '',
	# defined for low complexity read
	frac_n			=> '', # fraction of base n/N
	frac_mono		=> '', # fraction of mono-nt
	frac_di			=> '', # fraction of any combination of 2 bases 	
);

GetOptions(
	"h"						=> \$option{h},
	"silent"			=> \$option{silent},
	# general settings 
	"p=i" 				=> \$option{p},  #number of processors to use
	"ibam=s"			=> \$option{ibam}, # input bam file
	"ifq=s" 			=> \$option{ifq},	 # input fastq pair 1
	"ifq2=s"			=> \$option{ifq2}, # input fastq pair 2
	"batch=i"			=> \$option{batch},
	"odir=s"			=> \$option{odir},
	"oprefix=s"		=> \$option{oprefix},
	# defined for low complexity read
	"frac_n=f"		=> \$option{frac_n}, # fraction of base n/N
	"frac_mono=f"	=> \$option{frac_mono}, # fraction of mono-nt
	"frac_di=f"		=> \$option{frac_di}, # fraction of any combination of 2 bases 		
) || printHelp (); 

if ($option{h}) { printHelp();}

my $ibam = $option{ibam};
my $ifq = $option{ifq};
my $ifq2 = $option{ifq2};
my $batch = $option{batch};
my $odir = $option{odir};
my $oprefix = $option{oprefix};
my $silent_option = "";
if ($option{silent}) { $silent_option = " -silent "; }

unless ($ibam || ($ifq && $ifq2)) { print "\n\tERR: either -ibam or -ifq and -ifq2 should be specified\n"; printHelp ();}
unless ($odir && $oprefix) {	print "\n\tERR: -odir and -oprefix should be specified\n"; printHelp ();}

 
if ($ibam) {
#	$ibam =~ /.+\/(.+?)\..+/;
#	$prefix = $1;
	$ifq = $odir."/$oprefix.1.fq";
	$ifq2 = $odir."/$oprefix.2.fq";
	#print $prefix."\n".$ifq."\n".$ifq2."\n";
	
	# convert the BAM file to fastq format 
	bam2fastq ($ibam, $ifq, $ifq2);
} 

# ------------------------- convert fastq to sortable format ----------------------------------------------
#my $debug = 1;
#if ($debug != 1) {

my $begin_time = Time::HiRes::gettimeofday();
my $start_time = Time::HiRes::gettimeofday();

my $otsv = $odir."/$oprefix.tsv";

unless ($option{silent}) { print "\n\tConvert\n\t$ifq \n\tand \n\t$ifq2 \n\tto $otsv\n\n"; }

open (FQ, $ifq) or die $!;
open (FQ2, $ifq2) or die $!;
open my $fh_tsv, '>', $otsv or die "...$!";
my $readpair_cnt = 0;
my $idx = 0;
my $batch_idx = 0;
my @batch_records; 
my @record;					# to go from here...
while ((my $line = <FQ>) && (my $line2 = <FQ2>)) {
	chomp ($line);
	chomp ($line2);	
	if ($idx == 0) {
		if ($line =~ /^\@/ && $line2 =~ /^\@/) {
			++ $readpair_cnt;
			$idx = 1;
			my @entries = split (/\s/, $line);
			my $header = $entries[0];
			chomp ($header);
			push (@record, $header);
		} else { die "potential error in $ifq or $ifq2 format\n";	}	
	} elsif ($idx == 1) {
		push (@record, $line.$line2);
		$idx = 2;
	} elsif ($idx == 2) {
		$idx = 3;
	} elsif ($idx == 3) {
		push (@record, $line.$line2);		
		push (@batch_records, [ @record ]);
		@record = ();
		$idx = 0;
		++ $batch_idx;
		if ($batch_idx >= $batch) {
			print_2darray ($fh_tsv, \@batch_records);
			@batch_records = ();
			$batch_idx = 0;
		}
	}
}

if ($batch_idx > 0) {
		print_2darray ($fh_tsv, \@batch_records);	
		@batch_records = ();
		$batch_idx = 0;		
} 
close (FQ);
close (FQ2);
close ($fh_tsv);

unless ($option{silent}) { print "\t\tnumber of total read pairs: $readpair_cnt\n"; }

my $end_time = Time::HiRes::gettimeofday();

unless ($option{silent}) { printf("\n\t\tTime used: %.2f min\n", ($end_time - $start_time)/60); }
$start_time = $end_time;

#} # if (debug != 1)

# ------------------------- unix sort ----------------------------------------------
# first clear out all previously temporary stored sorted files
opendir (DIR, $odir) or die $!;	
while (my $file = readdir(DIR)) {
		if ($file =~ /_$oprefix/) { 
			my $path =  $odir."/$file";
			system ("rm $path"); 
		}
}
closedir DIR;

my $merged = $odir."/$oprefix.sorted.tsv";
unless ($option{silent}) { print "\n\tSplit, sort and merge $otsv to $merged\n\n"; }

#if ($debug != 1) {

#my $readpair_cnt; # debug to remove

my $num_entry_per_file = $readpair_cnt/$option{p};
$num_entry_per_file = ($num_entry_per_file == int ($num_entry_per_file)) ? 
							$num_entry_per_file : int ($num_entry_per_file + 1); # round up
my $splitfile_prefix = $odir."/_"."$oprefix";

#my $otsv; # debug to remove

unless ($option{silent}) { print "\t\tSplit $otsv to smaller files with about $num_entry_per_file entires each\n\n"; }
system ("split -l $num_entry_per_file -a 3 $otsv $splitfile_prefix") == 0 or die "Split failed\n";

my @files_to_sort;
opendir (DIR, $odir) or die $!;	
# first get all files 
while (my $file = readdir(DIR)) {
		if ($file =~ /_$oprefix/) { push (@files_to_sort, $file); }
}
closedir DIR;
my @sorted_files;
foreach my $file (@files_to_sort) {
	my $path = $odir."/$file";
	my $opath = $odir."/$file.sorted";
	push (@sorted_files, $opath);
	unless ($option{silent}) { print "\t\tsort +1 -2 -T $odir -S 1G $path -o $opath &\n"; }
	system ("sort +1 -2 -T $odir -S 1G $path -o $opath &");
}


# ----------- wait till sorting of files are completed --------------------------------
foreach my $file (@sorted_files) {
	while (1) {
		if (-e $file) {	last;	} 
		else { sleep 5;	}
	}	
}

# clean files to sort
foreach my $file (@files_to_sort) {
	my $path = $odir."/$file";
	system ("rm $path");
}

#} # if ($debug != 1) 

# ------------------------ merge sorted files --------------------------------

my $common_str = $odir."/_$oprefix*.sorted";
unless ($option{silent}) { print "\n\t\tsort +1 -2 -m $common_str -o $merged\n"; } # debug to uncommon 
system ("sort +1 -2 -m $common_str -o $merged"); # debug to uncommon

$end_time = Time::HiRes::gettimeofday();
unless ($option{silent}) { printf("\n\t\tTime used: %.2f min\n", ($end_time - $start_time)/60); }
$start_time = $end_time;

# clean individual sorted files
foreach my $file (@sorted_files) {
	system ("rm $file");
}

# ----------- dupl removal, low complexity read pair removal --------------------------------
# method: for low complexity read pair, if one read survive, directly push into batch_records
# 				if both survive, push into batch_records, and retain it as @pre_entry to be compared
#					with the later read pairs

my $fq_output = $odir."/".$oprefix.".clean.1.fq";
my $fq2_output = $odir."/".$oprefix.".clean.2.fq";
my $singleton_output = $odir."/".$oprefix.".clean.s.fq";
my $start_time = Time::HiRes::gettimeofday();

my $frac_option;
if ($option{frac_n}) {
	$frac_option = $frac_option." -freq_n $option{frac_n} ";
}
if ($option{frac_mono}) {
	$frac_option = $frac_option." -freq_mono $option{frac_mono} ";
}
if ($option{frac_di}) {
	$frac_option = $frac_option." -freq_di $option{frac_di} ";
}

unless ($option{silent}) { 
	print "\n\t\t$duprm_c -itsv $merged -ofq $fq_output -ofq2 $fq2_output -os $singleton_output $frac_option\n\n"; 
}
system ("$duprm_c -itsv $merged -ofq $fq_output -ofq2 $fq2_output -os $singleton_output $frac_option $silent_option");

my $end_time = Time::HiRes::gettimeofday();
unless ($option{silent}) { printf("\n\t\tTime used: %.2f min\n\n", ($end_time - $start_time)/60); }
 
# ----------- clean: remove tmp files -------------------
print "\tclean up ...\n";

system ("rm $otsv");
system ("rm $merged");

if ($ibam) {
	system ("rm $ifq");
	system ("rm $ifq2");
}

unless ($option{silent}) { printf("\n\tTotal time for duprm: %.2f min\n\n", ($end_time - $begin_time)/60); }

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------- sub-functions --------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
 
###############################################################

sub printCmd {
	print "\nRunning cmd:\n";
	print "perl $0 -p $option{p}";
	if ($ibam) {	print " -ibam $ibam"; } 
	else {
		print " -ifq $ifq -ifq2 $ifq2";
	}
	print " -odir $odir";
	print "\n\n";
}
###############################################################

sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Given a BAM or paired-fastq files as input, concatenate paired reads into a TSV format file. Sort this file then\n";
		print "\tremove low complexity and duplicated read-pairs.\n\n";
		print "usage: ./$0 {-ibam [input.bam] or -ifq [fwd.fq] -ifq2 [rv.fq]} -odir [scratch_folder_for_sort] -oprefix [output_prefix]\n\n";
		print "-silent: opt; no print out to screen\n";
		print "-p: default 8; number of threads to be used\n";
		print "-ibam: input BAM file; needs to be specified if -ifq and -ifq2 were not\n";
		print "-ifq & -ifq2: input fastq PE files; needs to be specified if -ibam is not\n";
		print "-batch: default 0.5mil; number of read pairs to process\n";
		print "-odir: scratch space for sorting\n";		
		print "-oprefix: output prefix name\n\n";
		print "Parameters for specifying low complexity read\n";
		print "-frac_n: default 0.4; fraction of base n/N\n";
		print "-frac_mono: default 0.5; fraction of mono-nt\n";
		print "-frac_di: default 0.8; fraction of any combination of 2 bases\n\n";				
		print "NOTE: require Perl 5.10+; Assuming reads consists of either {a,c,g,t,n} or {A, C, G, T, N} \n\t
		the output will be named as $oprefix.1.clean.fq and $oprefix.2.clean.fq\n\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

###############################################################
#brief	check if a read is low complexity
sub is_low_complexity {
	my $read = shift;
	my $counter_ref = shift;
	my $frac_ref = shift;

	my $len = length ($read);
	for (my $i = 0; $i < $len; ++ $i) {
		my $base = substr ($read, $i, 1);
		uc($base);
		if ($base eq 'A') {++ $counter_ref->[0];} 
		elsif ($base eq 'C') { ++ $counter_ref->[1]; }
		elsif ($base eq 'G') { ++ $counter_ref->[2]; }
		elsif ($base eq 'T') { ++ $counter_ref->[3]; }
		elsif ($base eq 'N') { ++ $counter_ref->[4]; }
		else { print "Unknown char: $base found in read $read\n"; }
	}
	
	my @sorted_cnts = sort {$b <=> $a} @$counter_ref;
	my $n_frac = $counter_ref->[4]/$len;
	my $mono_frac = $sorted_cnts[0]/$len;
	my $di_frac = ($sorted_cnts[0] + $sorted_cnts[1])/$len;
	
#	print $n_frac."\t".$mono_frac."\t".$di_frac."\n";
	if ($n_frac >= $frac_ref->[0] || $mono_frac >= $frac_ref->[1] 
		|| $di_frac>= $frac_ref->[2]) {
		return 1;		
	}
	return 0;
} # is_low_complexity

###############################################################
# print @records = (header, read, read2, qual, qual2) tsv to paired fastq files
sub print_2fq {
	my $fq_hd = shift;
	my $fq2_hd = shift;
	my $records_ref = shift;
	
	foreach my $oned (@$records_ref){
		print $fq_hd $$oned[0]."/1\n";
		print $fq_hd $$oned[1]."\n";
		print $fq_hd "+\n";
		print $fq_hd $$oned[3]."\n";

		print $fq2_hd $$oned[0]."/2\n";
		print $fq2_hd $$oned[2]."\n";
		print $fq2_hd "+\n";
		print $fq2_hd $$oned[4]."\n";		
	}
}
###############################################################

# print 2d array to file
sub print_2darray {
	my $fh = shift;
	my $tdarray_ref = shift;
	foreach my $oned (@$tdarray_ref){
		foreach my $twod (@{$oned}) {
			print $fh $twod."\t";
		}
		print $fh "\n";
	}
} # print_2darray 

###############################################################
# convert each the BAM file to fastq format, the resulting file contains specified label $label 
sub bam2fastq {
	my $bam_file = shift;
	my $ofq = shift;
	my $ofq2 = shift;
	
	# set NON_PF=True to use all reads instead of just PF reads
		
	print "\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=SILENT INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2\n\n";
	system ("\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=SILENT INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2");  	
} # bam2fastq


#-----------------------------------------------------------------------------------
# perl version of duplicate removal -- this has been replaced by c++ version
# ----------------------------------------------------------------------------------
#open my $fq_fh, '>', $fq_output or die "...$!";
#open my $fq2_fh, '>', $fq2_output or die "...$!";
#open (SORTEDTSV, $merged) or die $!;
#my $cnt_surviving = 0;
#my $cnt_low_complexity_read = 0;
#
#my @pre_counter = (0,0,0,0,0); # read 1 {acgtn/ACGTN}
#my @pre_counter2 = (0,0,0,0,0); # read 2 {acgtn/ACGTN}
#my @pre_entry = (); # (header, read, read2, qual, qual2)
#my $debug_counter = 0; 
#
#while (my $line = <SORTEDTSV>) {
#	++ $debug_counter; 
#	if ($debug_counter % 10000 == 0) {
#		print "counter = $debug_counter\t$cnt_surviving\t$cnt_low_complexity_read\n";
#	}
#	my @cur_counter = (0,0,0,0,0);
#	my @cur_counter2 = (0,0,0,0,0);
#	
#	my @entries = split ('\t', $line); 
#	my $len = length($entries[1]);
#	my $read = substr ($entries[1], 0, $len/2);
#	my $read2 = substr ($entries[1], $len/2, $len/2);
# 
#	my $is_lc = is_low_complexity ($read, \@cur_counter, \@frac);
#	my $is_lc2 = is_low_complexity ($read2, \@cur_counter2, \@frac);
#	
#	if ($is_lc == 1) { ++ $cnt_low_complexity_read; }
#	if ($is_lc2 == 1) { ++ $cnt_low_complexity_read; }
#
#	unless ($is_lc == 1 && $is_lc2 == 1) { # not both are low complexity
#			if ($is_lc == 1) { # first read is low complexity
#					print $fq_fh $entries[0]."/1\nN\n+\n#\n";
#					print $fq2_fh $entries[0]."/2\n$read2\n+\n".substr ($entries[2], $len/2, $len/2)."\n";
#					++ $cnt_surviving;
#				
#			} elsif ($is_lc2 == 1) { # 2nd read is low complexity
#					print $fq_fh $entries[0]."/1\n$read\n+\n".substr ($entries[2], 0, $len/2)."\n";
#					print $fq2_fh $entries[0]."/2\nN\n+\n#\n";
#					++ $cnt_surviving;	
#			} else { # check if it is the same compared to the previous entry
#
#				unless (@pre_entry) { # first valid read pair; directly register as @pre_entry and push into batch_records		
#					@pre_entry = ($entries[0], $read, $read2, substr ($entries[2], 0, $len/2),
#												 substr ($entries[2], $len/2, $len/2));
#					@pre_counter = @cur_counter;
#					@pre_counter2 = @cur_counter2;							 
#					print $fq_fh $entries[0]."/1\n$read\n+\n".substr ($entries[2], 0, $len/2)."\n";
#					print $fq2_fh $entries[0]."/2\n$read2\n+\n".substr ($entries[2], $len/2, $len/2)."\n";
#					$cnt_surviving += 2;
#				} else { # check duplicity by comparing with the @pre_entry 	
#					if ((@cur_counter ~~ @pre_counter) && (@cur_counter2 ~~ @pre_counter2)) { # compare nt counts first for speeding up
							# use string comparison to validate 
#						unless (($read eq $pre_entry[1]) && ($read2 eq $pre_entry[2])) { # not equal 
#							@pre_entry = ($entries[0], $read, $read2, substr ($entries[2], 0, $len/2),
#														 substr ($entries[2], $len/2, $len/2));
#							@pre_counter = @cur_counter;
#							@pre_counter2 = @cur_counter2;							 														 

#							print $fq_fh $entries[0]."/1\n$read\n+\n".substr ($entries[2], 0, $len/2)."\n";
#							print $fq2_fh $entries[0]."/2\n$read2\n+\n".substr ($entries[2], $len/2, $len/2)."\n";

#							$cnt_surviving += 2;				
#						}
												
#					}	else { # not equal 

#						@pre_entry = ($entries[0], $read, $read2, substr ($entries[2], 0, $len/2),
#													 substr ($entries[2], $len/2, $len/2));
#						@pre_counter = @cur_counter;
#						@pre_counter2 = @cur_counter2;							 
													 
#						print $fq_fh $entries[0]."/1\n$read\n+\n".substr ($entries[2], 0, $len/2)."\n";
#						print $fq2_fh $entries[0]."/2\n$read2\n+\n".substr ($entries[2], $len/2, $len/2)."\n";
#						$cnt_surviving += 2;
#					}			
#				}
#			}
#	} # unless ($is_lc == 1 && $is_lc2 == 1) 

#}

#close (SORTEDTSV);
#close ($fq_fh);
#close ($fq2_fh);
