# Given read names, specified in a tsv file as the 1st field (e.g. bowtie output),
# remove any matching reads in the input fastq paired or single end files. output
# as fastq paired files (if both end remain) + a single end fastq file containing
# only single end reads

#!/usr/bin/perl
use strict;
use Getopt::Long;

my %option = (
	h         => '',
	silent 		=> '',
	# I/O setting 
	itsv			=> '', 
	ipfqlist	=> '',
	isfqlist	=> '',
	opfqlist	=> '',
	osfq			=> '',
);

GetOptions(
	"h"							=> \$option{h},
	"silent"				=> \$option{silent},
	"itsv=s"				=> \$option{itsv},
	"ipfqlist=s"		=> \$option{ipfqlist},
	"isfqlist=s" 		=> \$option{isfqlist},
	"opfqlist=s"		=> \$option{opfqlist},
	"osfq=s" 				=> \$option{osfq},
) || printHelp (); 

if ($option{h}) { printHelp();}

my $itsv = $option{itsv};
my $ipfqlist = $option{ipfqlist};
my $isfqlist = $option{isfqlist};
my $opfqlist = $option{opfqlist};
my $osfq = $option{osfq};
my $max_frag_len = $option{max_frag_len};

unless ($itsv && $osfq) {
	print "\nErr: both -itsv and -osfq should be specified\n";
	printHelp ();
}

if (!$ipfqlist && !$isfqlist) {
	print "\nErr: at least one of -ipfqlist and -isfqlist should be specified\n";
	printHelp ();
}

if (($ipfqlist && !$opfqlist) || (!$ipfqlist && $opfqlist)) {
	print "\nErr: both or none of -ipfqlist and -opfqlist should be specified\n";
	printHelp ();
}

my @ipfqfiles = split (',', $ipfqlist);
my $num_ipfq = $#ipfqfiles + 1;
my @isfqfiles = split (',', $isfqlist);
my $num_isfq = $#isfqfiles + 1;
my @opfqfiles = split (',', $opfqlist);
my $num_opfq = $#opfqfiles + 1;

if ($num_ipfq) {
	if ($num_ipfq % 2 != 0) {
		print "\nErr: odd number $num_ipfq ipfqlist files found\n";
		printHelp ();
	}
}

if ($num_opfq != $num_ipfq) {
		print "\nErr: $num_opfq files in opfqlist != $num_ipfq files in ipfqlist\n";
		printHelp ();
}

my $min_len = 1000000000;
my $max_len = 0;
my $num_reads = 0;
my $sum_len;

# ------------- read through $itsv file -------------
unless ($option{silent}) {
	print "\n\tAnalyze itsv file: $itsv\n\n";
}
my %rname_hash;
open (ITSV, $itsv) or die $!;
while (my $line = <ITSV>) {
	next if ($line =~ /^\s*$/); #skip blank lines
	chomp($line);
	my @entries = split ('\t', $line);
	$rname_hash{$entries[0]} = 1;	
	++ $num_reads;
}
close (ITSV);
print "\tnum_reads in $itsv: $num_reads\n\n";

# -----------------------------------------------------------

open my $osfq_fh, ">$osfq" or die $!; # output single fq file

# ---------- search in paired fastq files ------------------
for (my $i = 0; $i < $num_ipfq; $i += 2) {

	unless ($option{silent}) {
		print "\tParse paired-end files: $ipfqfiles[$i] and $ipfqfiles[$i+1]\n\n";
	}
	
	open (FQ, $ipfqfiles[$i]) or die $!;
	open (FQ2, $ipfqfiles[$i+1]) or die $!;
	open my $ofq_fh, ">$opfqfiles[$i]" or die $!;
	open my $ofq2_fh, ">$opfqfiles[$i+1]" or die $!;

	while ((my $line = <FQ>) && (my $line2 = <FQ2>)) {
		my @first_end;
		my @second_end;
		
		if ($line =~ /^\@/ && $line2 =~ /^\@/) { 
				push (@first_end, $line);
				push (@second_end, $line2);
				
				for (my $i = 0; $i < 3; ++ $i) {
						$line = <FQ>;
						$line2 = <FQ2>;
						push (@first_end, $line);
						push (@second_end, $line2);						
				}	
				
#				print "@first_end\n@second_end\n";

				my $header = $first_end[0];
				my $header2 = $second_end[0];
				chomp ($header);
				chomp ($header2);
				$header = substr($header, 1, length($header) - 1); # remove the first char of readname to be consistent with bowtie output 
				$header2 = substr($header2, 1, length($header2) - 1); 
				
				# just to comment out the following two lines if input fastq pairs are not aligned by bowtie in paired-end fashion
#				$header = $header."/1";  # adding /1 and /2 to be consistent with bowtie paired hits
#				$header2 = $header2."/2";

				#	new strategy: if one end if found, remove both
				unless (exists $rname_hash{$header} || exists $rname_hash{$header2}) { # first or second end found
						output_read (\@first_end, $ofq_fh);
						output_read (\@second_end, $ofq2_fh);									
				}	
								
#				if (exists $rname_hash{$header}) { # first end found
#					unless (exists $rname_hash{$header2}) { # second end not found
#						output_read (\@second_end, $osfq_fh);
#					}
#				} else {
#					if (exists $rname_hash{$header2}) { # second end found
#						output_read (\@first_end, $osfq_fh);
#					} else {
#						output_read (\@first_end, $ofq_fh);
#						output_read (\@second_end, $ofq2_fh);
#					}					
#				}

		} else { 
				die "\npotential error in the format of file $ipfqfiles[$i] or $ipfqfiles[$i+1]\n\n";	
		}	 
	} # while 
	
	close (FQ);
	close (FQ2);
	close (ofq_fh);
	close (ofq2_fh);
} # for 


# ---------- search in single-end fastq files ------------------

for (my $i = 0; $i < $num_isfq; $i += 2) {

	unless ($option{silent}) {
		print "\tParse single end file: $isfqfiles[$i]\n\n";		
	}
	open (FQ, $isfqfiles[$i]) or die "Cannot open $isfqfiles[$i]\n";
	while (my $line = <FQ>) {
		my @first_end;
		if ($line =~ /^\@/) {
			push (@first_end, $line);
			for (my $i = 0; $i < 3; ++ $i) {
				$line = <FQ>;
				push (@first_end, $line);
			}
			my $header = $first_end[0];
			chomp ($header);
			$header = substr($header, 1, length($header) - 1); # remove the first char of readname to be consistent with bowtie output 
			unless (exists $rname_hash{$header}) { # not found
				output_read (\@first_end, $osfq_fh);
			}	
		} else {
			die "\npotential error in the format of file $ipfqfiles[$i]\n\n";	
		}
	}
	close (FQ);
}	
close (osfq_fh);
################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Given read names specified in a tsv file as the 1st field (e.g. bowtie output),\n";
		print "\tremove any matching reads in the input fastq paired or single end files. Output\n";
		print "\tas fastq paired files (if both end remain) + a single end fastq file containing\n";
		print "\tall single end reads\n\n";
		print "usage: ./$0 -itsv [names.tsv] -ipfqlist [1.fwd.fq],[1.rv.fq]... -isfqlist [a.fq],[b.fq]... -opfqlist [1.fwd.clean.fq][1.rv.clean.fq]... -osfq [clean.fq]... \n\n";
		print "-h: help; for cmd options\n";
		print "-silent: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-itsv: input file in tsv format with 1st field being read names\n";
		print "-ipfqlist: comma separated even number of input fastq file list, where adjacent two are paired-ends\n";
		print "-isfqlist: input single end fastq file list\n";
		print "-opfqlist: output fastq file list, where adjacent two are paired-ends (same number as -ipfqlist)\n";
		print "-osfq: output single end fastq file\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}

# output fastq read entry to a fastq file handle
sub output_read {
	my $array_ref = shift;
	my $fh = shift;
	foreach (@$array_ref) {
		print $fh $_;
	}
}