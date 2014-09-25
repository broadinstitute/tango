# @description: select from SAM file alignment conforming criteria:
#	i) alignment/input_len >= min_perc_span, and ii) similarity >= min_perc_sim;

#!/usr/bin/perl

use strict;
use Getopt::Long;

my $bwa_cmd = "bwa mem -t 8 -k 23 -w 4 -T 70 /seq/viral/analysis/xyang/FUO/DB/bwa_meta /seq/viral/analysis/xyang/programs/M-Vicuna/input/338.2m.fwd.fq";
my $min_perc_span = 0.95; #$ARGV[1];
my $min_perc_sim = 0.95; #$ARGV[2];
open(CMD, "$bwa_cmd |") or die $!;

while (my $line = <CMD>)
{

	my @entries = split ('\t', $line);
	if ($entries[0] =~ /^@/) {		
		# print $line."\n";
	} else {
 
		if ($#entries + 1 < 11) {
			print "[error] line: $line has < 11 entries\n";
			exit;
		}
		
		my $query_name = $entries[0];
		my $hit_name = $entries[2];
		my $hit_pos = $entries[3];
		my $cigar = $entries[5];
		my $query_seq = $entries[9];
		
		#$cigar = "10P2M1I2D45M";
		
		my $num_match = get_cigar_char_cnt ($cigar, 'M');
		next if $num_match == 0;
		my $query_len = length ($query_seq);
		
		my $query_seq_len = $query_len + get_cigar_char_cnt ($cigar, 'H');	# raw sequence input length
		my $span_len = $query_len - get_cigar_char_cnt ($cigar, 'S'); # minus soft-clip
		
		my $perc_span = $span_len * 100 / $query_seq_len;
		my $perc_sim =  $num_match * 100 / get_cigar_char_cnt ($cigar, "[M,I,D,N,=,X]");
		
		#{ # debug code
		#	print $query_seq."\tlen=$query_len\n".$cigar."\t".
		#	$num_match."\t".$num_span_len."\t".$num_hard_clip."\n".$perc_span."\t".$perc_sim."\n";
		#	exit;
		#}
		
		if (($perc_span >= $min_perc_span) && ($perc_sim >= $min_perc_sim)) {
			print $query_name."\t".$hit_name."\t"."$perc_span\t$perc_sim\n";		
 		}		
	} #
	
}

close(CMD);

###############################################################

# return the number of specific cigar character in a cigar string 
sub get_cigar_char_cnt {
	my $cigar = shift;
	my $c = shift;
	
	my $total = 0;
	my @entries = split ($c, $cigar);
	foreach my $elem (@entries) {
		$elem =~ /(\d+)$/;
		$total += $1;
	}	
	return $total;
} # get_cigar_char_cnt