#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use warnings;

my %option = (  
	h				=> '',
	quiet			=> '',
	p				=> 10,
	# --------  input fastq pairs ---------
	fq				=> '',
	fq2			=> '',
	# -------- aligner ---------------
	blast			=> '',
	bt2			=> '',			
	# ------- input ---------------
	fa				=> '',
  	odir		    => '',
	oprefix		=> '',
);

GetOptions(
	"h"		    => \$option{h},
	"quiet"		=> \$option{quiet},
	"p=i" 		=> \$option{p},
	# -------- aligner ---------------
	"blast"		=> \$option{blast},
	"bt2"		=> \$option{bt2},
	# -----------input---------
  "fa=s"			=> \$option{fa},
  "fq=s"			=> \$option{fq},
  "fq2=s"		=> \$option{fq2},
  
  "odir=s"      => \$option{odir},
  "oprefix=s"	=> \$option{oprefix},
) || die("Problem processing command-line options: $!\n");

my $oprefix = $option{oprefix};
my $odir = $option{odir};
unless ($oprefix && $odir) {
	print "-oprefix and -odir need to be specified\n";
	printHelp();
}

my $is_bt2 = $option{bt2};
my $is_blast = $option{blast};

if (($is_bt2 && $is_blast) || (! $is_bt2 && ! $is_blast)) {
	print "one and only one of -is_bt2 and -is_blast should be specified\n";
	printHelp();
}

my $metaphlan_py = "/seq/viral/analysis/xyang/external_programs/metaphlan1.7.8/metaphlan.py";
my $bowtie2 = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/bowtie2_2.1.0/bowtie2";
my $bowtie2db = "/seq/viral/analysis/xyang/external_programs/metaphlan1.7.8/bowtie2db/mpa";
my $blast = "/xchip/koch/softwares/ncbi-blast-2.2.27+/bin/blastn";
my $blastdb = "/seq/viral/analysis/xyang/external_programs/metaphlan1.7.8/blastdb/mpa";
my $fq2fa_pl = "/seq/viral/analysis/xyang/scripts/others/fastq2fasta.pl";


# read input of fq or fa 
my $fa = $option{fa};
my $fq = $odir."/$oprefix.fq";
my $output = $odir."/$oprefix.phy.txt";

unless ($fa) { 
	$fa = $odir."/$oprefix.fa"; 
	if (-e $fa) { system ("rm $fa"); }
}

if ($is_blast) { # blast as aligner 
	if ($option{fq}) {
		my $tmp_fa = $odir."/$oprefix.tmp";
	
		print "	[CMD] perl $fq2fa_pl $option{fq} $tmp_fa\n";

		system ("perl $fq2fa_pl $option{fq} $tmp_fa");
		print "	[CMD] cat $tmp_fa.fa >> $fa\n";

		system ("cat $tmp_fa.fa >> $fa"); 		
		if ($option{fq2}) {
			print "	[CMD] perl $fq2fa_pl $option{fq2} $tmp_fa\n";
			system ("perl $fq2fa_pl $option{fq2} $tmp_fa");
			print "	[CMD] cat $tmp_fa.fa >> $fa\n";
			system ("cat $tmp_fa.fa >> $fa"); 		
		}
	
		system ("rm $tmp_fa.fa");
		system ("rm $tmp_fa.qual");
	}
	my $blast_out = $odir."/$oprefix.blast.txt";
	print "	[CMD] $metaphlan_py --word_size 12 --blastn_exe $blast --blastdb $blastdb --nproc $option{p} --blastout $blast_out $fa $output\n";
	system ("$metaphlan_py --word_size 12 --blastn_exe $blast --blastdb $blastdb --nproc $option{p} --blastout $blast_out $fa $output");
	
} else { # bowtie2 
	if (-e $fq) { 
		print "\n\t[CMD] rm $fq\n";		
		system ("rm $fq"); 
	}
	if ($option{fq}) { 
		print "\n\t[CMD] cat $option{fq} >> $fq\n";
		system ("cat $option{fq} >> $fq"); 
	}
	if ($option{fq2}) { 		
		print "\n\t[CMD] cat $option{fq2} >> $fq\n";
		system ("cat $option{fq2} >> $fq"); 
	}
	
	my $bt_out = $odir."/$oprefix.bt2.txt";
	if (-e $bt_out) { system ("rm $bt_out"); }
	
	print "\n\t[CMD] $metaphlan_py --bowtie2_exe $bowtie2 --bowtie2db $bowtie2db --bt2_ps sensitive-local --nproc $option{p} --bowtie2out $bt_out $fq $output\n";
	system ("$metaphlan_py --bowtie2_exe $bowtie2 --bowtie2db $bowtie2db --bt2_ps sensitive-local --nproc $option{p} --bowtie2out $bt_out $fq $output");
	
} 

#/xchip/koch/softwares/metaphlan/metaphlan.py --min_cu_len 10000 < Step 1 output file file extension .outfmt6.txt > sample.profiling_output.txt


################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Run metaphlan\n";
		print "\nusage: ./$0 -ifq [fwd.fq] -ifq2 [rv.fq] -oprefix [oprefix]\n\n";
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n";
		print "-p: default 10, number of proc\n";
		print "-bt2: use bowtie2 as aligner, mutually exclusive with -blast\n";
		print "-blast: use blastn as aligner, mutulally exclusive with -bt2\n";
		print "-fa: input fasta file\n";
		print "-fq: input fq file\n";
		print "-fq2: second fq file\n";
		print "-odir: output dir\n";
		print "-oprefix: output prefix\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}
