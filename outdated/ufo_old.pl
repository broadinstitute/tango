# @brief
# Prerequisite: 
# 1) Perl 5.10 +
# 2) set script path in ufo.pl
# 3) set script path in duprm_by_sort.pl, 
# 4) gcc-4.7.2 or above

#!/usr/bin/perl

use strict;
use Getopt::Long;

	#--------------------------------------------------------------------------------------------------
	# set dependent programs/scripts/data path
	# system 
	my $perl_pl = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/perl_5.10.1/bin/perl";

	# dependent programs 
	my $sam2fastq_jar = "/seq/software/picard/current/bin/SamToFastq.jar";
	my $trimmomatic_jar = "/seq/viral/analysis/xyang/external_programs/Trimmomatic-0.29/trimmomatic-0.29.jar";
	my $trimref = "/seq/viral/analysis/xyang/FUO/DB/TruSeq3-PE.fa"; # Trimmomatic dependent
	my $bowtie_c = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/bowtie_0.12.9/bowtie";

	# our package
	#my $trimref = "/seq/viral/analysis/xyang/FUO/DB/adaptor_primer.fasta";
	#my $duprm_pl = "/seq/viral/analysis/xyang/FUO/scripts/duprm_by_sort.pl";
	my $duprm_c = "/seq/viral/analysis/xyang/programs/DupRm/bin/duprm"; 
	my $rmhits_pl = "/seq/viral/analysis/xyang/FUO/scripts/rm_hits.pl";
	my $mvicuna_c = "/seq/viral/analysis/xyang/programs/M-Vicuna/bin/mvicuna";

	# reference database
	my $db_hrna = "/seq/viral/analysis/xyang/FUO/DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA";
	my $db_hg = "/seq/viral/analysis/xyang/FUO/DB/bowtie_idx_hg19";
	my $db_plas = "/seq/viral/analysis/xyang/FUO/DB/bowtie_idx_plasmodium";
	my $db_meta = "/seq/viral/analysis/xyang/FUO/DB/bowtie_idx_metagenomics_contaminants_v2";
	#--------------------------------------------------------------------------------------------------

	my %option = (
		h           		=> '',
		silent 				=> '',
		# I/O setting 
		ibams				=> '', # input comma separated BAM files, where each denotes a paired fq files
		ipfqs 				=> '',	 # input paired comma separated fastq files 
		odir				=> '',  
		# performance setting 
		p 					=> 8,		 # performance control; num of cores to use
		batch				=> 500000,   # performance control; batch of reads/read-pairs to be loaded in memory

		# low complexity read setting 
		frac_n			=> '', # fraction of base n/N
		frac_mono		=> '', # fraction of mono-nt
		frac_di			=> '', # raction of any combination of 2 bases

		# Trimmomatic setting
	#	trim_mqual  => 25,
	#	trim_w			=> 4,	# trimming window, average qual needs to be >= $trim_mqual 
	#	trim_mlen		=> 70,	# min read length to keep post-trimming
		trim_mqual  => 15,
		trim_w			=> 10,	# trimming window, average qual needs to be >= $trim_mqual 
		trim_mlen		=> 50,	# min read length to keep post-trimming
	
		#bowtie setting
		bt_max_frag_len => '',
	);

	GetOptions(
		"h"					=> \$option{h},
		"silent"				=> \$option{silent},
		# I/O setting
		"ibams=s"			=> \$option{ibams}, 	
		"ipfqs=s" 			=> \$option{ipfqs},	
		"odir=s"			=> \$option{odir}, 
	
		# performance setting
		"p=i" 				=> \$option{p},  #number of processors to use
		"batch=i"			=> \$option{batch},   

		# low complexity read setting
		"frac_n=f"		=> \$option{frac_n}, # fraction of base n/N
		"frac_mono=f"	=> \$option{frac_mono}, # fraction of mono-nt
		"frac_di=f"		=> \$option{frac_di}, # fraction of any combination of 2 bases 	
	
		# Trimmomatic setting
		"trim_mqual=i" 	=> \$option{trim_mqual},
		"trim_w=i"			=> \$option{trim_w},
		"trim_mlen=i"		=> \$option{trim_mlen},

		# bowtie setting 
		"bt_max_frag_len=i" => \$option{bt_max_frag_len},

	) || printHelp (); 

	if ($option{h}) { printHelp();}

	# ----------------------------------------------- get input parameters -----------------------------------------------------------
	my $ibams = $option{ibams};
	my $ipfqs = $option{ipfqs};
	my $odir = $option{odir};
	my $silent_option = "";

	if ($option{silent}) { $silent_option = " -silent "; }

	unless ($ibams || $ipfqs) { print "\n\tERR: either -ibams or -ipfqs should be specified\n"; printHelp ();}
	if ($ibams && $ipfqs) { print "\n\tERR: only one of -ibams and -ipfqs should be specified\n"; printHelp ();}
	unless ($odir) { print "\n\tERR:  -odir should be specified\n"; printHelp ();}

	# --------------------- convert each bam file to paired fastq files, and record in $ipfqs ----------------------------------------
	
	if ($ibams) { $ipfqs = generate_input_fastq ($ibams, $odir); } 
	
	# ---------------------------------------- prepare output paired fastq files -----------------------------------------------------
	my $drm_opfqs;
	my $prm_opfqs;
	my $prm_ofa;
	generate_mvicuna_output_filenames (\$drm_opfqs, \$prm_opfqs, \$prm_ofa, $ipfqs, $odir);
 
	# --------------------------------------------------------------------------------------------------------------------------------
	# -------------------------------- duplicate & low complexity read removal & paired read merging ---------------------------------
	# --------------------------------------------------------------------------------------------------------------------------------

	print "\n[[[De novo duplicate removal, low complexity read removal, paired-read merging]]]\n\n";
	my $frac_option;
	if ($option{frac_n}) {	$frac_option .= " -drm_max_n_perc $option{frac_n} "; }
	if ($option{frac_mono}) { $frac_option .= " -drm_max_mono_perc $option{frac_mono} "; }
	if ($option{frac_di}) { $frac_option .= " -drm_max_di_perc $option{frac_di} "; }

	print "\n$mvicuna_c -batch $option{batch} -pthreads $option{p} -pfq $ipfqs -t_drm -drm_ofq $drm_opfqs $frac_option -t_prm -prm_ofq $prm_opfqs -prm_ofa $prm_ofa\n\n"; 
	system ("$mvicuna_c -batch $option{batch} -pthreads $option{p} -pfq $ipfqs -t_drm -drm_ofq $drm_opfqs $frac_option -t_prm -prm_ofq $prm_opfqs -prm_ofa $prm_ofa"); 

	exit; 

my $oprefix;

	if ($ibams) {
		print "\nperl $duprm_pl -ibams $ibams -odir $odir -oprefix $oprefix -p $option{p} -batch $option{batch} $frac_option $silent_option\n\n";
		system("$perl_pl $duprm_pl -ibam $ibams -odir $odir -oprefix $oprefix -p $option{p} -batch $option{batch} $frac_option $silent_option") == 0 or die $!;

	} else {
		print "\nperl $duprm_pl -ipfqs $ipfqs -odir $odir -oprefix $oprefix -p $option{p} -batch $option{batch} $frac_option $silent_option\n\n";
		system("$perl_pl $duprm_pl -ipfqs $ipfqs -odir $odir -oprefix $oprefix -p $option{p} -batch $option{batch} $frac_option $silent_option") == 0 or die $!; 	
	}

	# output $odir."/".$oprefix.".clean.1.fq",  $odir."/".$oprefix.".clean.2.fq" and $odir."/".$oprefix.".clean.s.fq";

 
	# --------------------------------------------------------- Trimmomatic --------------------------------------------------------------------------
	print "\n[[[Trimming]]]\n\n";

	# Trimmomatic input 
	my $iclean_fq = $odir."/".$oprefix.".clean.1.fq";
	my $iclean_fq2 = $odir."/".$oprefix.".clean.2.fq";
	my $iclean_s = $odir."/".$oprefix.".clean.s.fq";
	unless ((-e $iclean_fq) && (-e $iclean_fq2) && (-e $iclean_s)) { print "\nERR: input files: 1) $iclean_fq or 2) $iclean_fq2 or 3) $iclean_s to trimmomatic doesn't exist \n"; exit; }

	# Trimmomatic output 
	my $trimfq = $odir."/".$oprefix.".trim.p1.fq";
	my $trimfq2 = $odir."/".$oprefix.".trim.p2.fq";
	my $trim_singleton = $odir."/".$oprefix.".trim.s1.fq";
	my $trim_singleton2 = $odir."/".$oprefix.".trim.s2.fq";
	my $trim_s = $odir."/".$oprefix.".trim.s.fq";

	# NOTE: putting the argument ILLUMINACLIP:$trimref:2:30:10 in front is more time consuming compared to when putting it at the end of
	# the cmd and the trimming results differ slightly; e.g. with 17mil paired-end reads, 50mins vs. 30mins for 8 cores on Hyacinth.
	# [-trimlog <logFile>]
	print "java -jar -Xmx4g $trimmomatic_jar PE -threads $option{p} -phred33 $iclean_fq $iclean_fq2 $trimfq $trim_singleton $trimfq2 $trim_singleton2 ILLUMINACLIP:$trimref:2:30:10 LEADING:$option{trim_mqual} TRAILING:$option{trim_mqual} SLIDINGWINDOW:$option{trim_w}:$option{trim_mqual} MINLEN:$option{trim_mlen}\n\n";
	system ("java -jar -Xmx4g $trimmomatic_jar PE -threads $option{p} -phred33 $iclean_fq $iclean_fq2 $trimfq $trim_singleton $trimfq2 $trim_singleton2 ILLUMINACLIP:$trimref:2:30:10 LEADING:$option{trim_mqual} TRAILING:$option{trim_mqual} SLIDINGWINDOW:$option{trim_w}:$option{trim_mqual} MINLEN:$option{trim_mlen}") == 0 or die $!;

	print "java -jar -Xmx4g $trimmomatic_jar SE -threads $option{p} -phred33  $iclean_s $trim_s ILLUMINACLIP:$trimref:2:30:10 LEADING:$option{trim_mqual} TRAILING:$option{trim_mqual} SLIDINGWINDOW:$option{trim_w}:$option{trim_mqual} MINLEN:$option{trim_mlen} \n\n";
	system ("java -jar -Xmx4g $trimmomatic_jar SE -threads $option{p} -phred33 $iclean_s $trim_s ILLUMINACLIP:$trimref:2:30:10 LEADING:$option{trim_mqual} TRAILING:$option{trim_mqual} SLIDINGWINDOW:$option{trim_w}:$option{trim_mqual} MINLEN:$option{trim_mlen}") == 0 or die $!;
	# clean trimmomatic input
	system ("rm $iclean_fq");
	system ("rm $iclean_fq2");
	system ("rm $iclean_s");

	# concatenation of all singletons
	print "cat $trim_singleton $trim_singleton2 >> $trim_s\n\n";
	system ("cat $trim_singleton $trim_singleton2 >> $trim_s");
	system ("rm $trim_singleton");
	system ("rm $trim_singleton2");

	# output $trim_s, $trimfq $trimfq2
	# ----------------------------------------- low complexity read removal post-trimming --------------------------------------------------------------
	print "\n[[[Low complexity read removal]]]\n\n";

	my $fq_output = $odir."/".$oprefix.".trim.lcrm.p1.fq";
	my $fq2_output = $odir."/".$oprefix.".trim.lcrm.p2.fq";
	my $s_output = $odir."/".$oprefix.".trim.lcrm.s.fq";
	my $s_tmp_output = $odir."/".$oprefix.".trim.lcrm.tmps.fq";
	# 1. deal w/ paired reads
	print "\n$duprm_c -ifq $trimfq -ifq2 $trimfq2 -ofq $fq_output -ofq2 $fq2_output -os $s_output $frac_option -rm_lc_only $silent_option\n\n"; 
	system ("$duprm_c -ifq $trimfq -ifq2 $trimfq2 -ofq $fq_output -ofq2 $fq2_output -os $s_output  $frac_option -rm_lc_only $silent_option");

	 # clean input
	system ("rm $trimfq"); 
	system ("rm $trimfq2");

	# 2. deal with singletons
	print "\n$duprm_c -ifq $trim_s -os $s_tmp_output  $frac_option -rm_lc_only $silent_option\n\n"; 
	system ("$duprm_c -ifq $trim_s -os $s_tmp_output  $frac_option -rm_lc_only $silent_option");

	# merge singletons
	print "cat $s_tmp_output >> $s_output\n\n";
	system ("cat $s_tmp_output >> $s_output");
	# clean input
	system ("rm $trim_s"); 
	system ("rm $s_tmp_output");

	# outputs: $fq_output, $fq2_output, $s_output

	# ----------------------------------------- Remove known reads --------------------------------------------------------------
	my $hit_hg = $odir."/".$oprefix.".hit.hg.txt";
	my $hit_hrna = $odir."/".$oprefix.".hit.hrna.txt";
	my $hit_meta = $odir."/".$oprefix.".hit.meta.txt";

	my $ipfqlist = $fq_output.",".$fq2_output;
	my $opfq = $fq_output.".clean.txt";
	my $opfq2 = $fq2_output.".clean.txt";
	my $opfqlist = $opfq.",".$opfq2;
	my $osfq = $s_output.".clean.txt";

	my $bowtie_prnt;
	if ($option{silent}) {
		$bowtie_prnt = " --quiet ";
	}


	#system ("$bowtie_c $db_meta -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -1 $fq_output -2 $fq2_output -X $option{bt_max_frag_len} $bowtie_prnt > $hit_meta");
	system ("$bowtie_c $db_meta -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq_output $bowtie_prnt > $hit_meta");
	system ("$bowtie_c $db_meta -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq2_output $bowtie_prnt >> $hit_meta");
	system ("$bowtie_c $db_meta -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $s_output $bowtie_prnt >> $hit_meta");
	system ("perl $rmhits_pl -itsv $hit_meta -ipfqlist $ipfqlist -isfqlist $s_output -opfqlist $opfqlist -osfq $osfq $silent_option");
	system ("mv $opfq $fq_output");
	system ("mv $opfq2 $fq2_output");
	system ("mv $osfq $s_output");

	#system ("$bowtie_c $db_hg -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -1 $fq_output -2 $fq2_output -X $option{bt_max_frag_len} $bowtie_prnt > $hit_hg");
	system ("$bowtie_c $db_hg -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq_output $bowtie_prnt > $hit_hg");
	system ("$bowtie_c $db_hg -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq2_output $bowtie_prnt >> $hit_hg");
	system ("$bowtie_c $db_hg -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $s_output $bowtie_prnt >> $hit_hg");
	system ("perl $rmhits_pl -itsv $hit_hg -ipfqlist $ipfqlist -isfqlist $s_output -opfqlist $opfqlist -osfq $osfq $silent_option");
	system ("mv $opfq $fq_output");
	system ("mv $opfq2 $fq2_output");
	system ("mv $osfq $s_output");

	#system ("$bowtie_c $db_hrna -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -1 $fq_output -2 $fq2_output -X $option{bt_max_frag_len} $bowtie_prnt > $hit_hrna");
	system ("$bowtie_c $db_hrna -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq_output $bowtie_prnt > $hit_hrna");
	system ("$bowtie_c $db_hrna -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $fq2_output $bowtie_prnt >> $hit_hrna");
	system ("$bowtie_c $db_hrna -p $option{p} -k 1 -v 2 --suppress 5,6,7,8 -q $s_output $bowtie_prnt >> $hit_hrna");
	system ("perl $rmhits_pl -itsv $hit_hrna -ipfqlist $ipfqlist -isfqlist $s_output -opfqlist $opfqlist -osfq $osfq $silent_option");
	system ("mv $opfq $fq_output");
	system ("mv $opfq2 $fq2_output");
	system ("mv $osfq $s_output");


 
	# ----------------------------------------- Remove Known Viral Reads --------------------------------------------------------------

	#################################################################################################################################
	# sub-functions 
	#############################################################################
	sub printCmd {
		print "Running cmd:\n";
		print "perl $0 -p $option{p} -ibams $ibams";
		print "\n";
	}

	#############################################################################
	sub printHelp {
			print "\n----------------------------------------------------------------------------\n";
			print "usage: perl $0 {-ibams [1.bam,2.bam,...] or -ipfqs [1.f.fq,1.r.fq,2.f.fq,2.r.fq...]} -odir [odir/] -oprefix [oprefix] \n\n";
		
			print "-h: help; for cmd options\n";
			print "-silent: no screen output for programs called\n\n";
		
			print "I/O setting\n";		
			print "-ibams: comma separated BAM files; mutually exclusive from -ipfqs \n";
			print "-ipfqs: comma separated paired fq files; mutually exclusive from -ibams\n";
			print "-odir: output directory\n";
			print "-oprefix: output file prefix\n\n";
		
			print "Performance setting\n";
			print "-p: default 8; number of threads\n";
			print "-batch: default 0.5 mil; number of reads/read-pairs to be in memory\n\n";
		
			print "duplicate and low complexity removal setting, used by mvicuna\n";
			print "-frac_n: default 0.4; fraction of base n/N\n";
			print "-frac_mono: default 0.5; fraction of mono-nt\n";
			print "-frac_di: default 0.8; fraction of any combination of 2 bases\n\n";
									
			print "Trimmomatic setting\n";
			print "-trim_mqual: default 20 (-phred33); min qual for prefix and suffix to be trimmed\n";
			print "-trim_w: default 4; any sliding window with avg min qual < trim_mqual will be trimmed\n";
			print "-trim_mlen: default 50; min read length to keep post-trimming\n";

			print "\n----------------------------------------------------------------------------\n";
			exit;
	}
	
	############################################################################
	# Given input list of BAM files, convert each to paired fastq files for each
	sub generate_input_fastq {
		my $ibams = shift;
		my $odir = shift;
		
		my $ipfqs;
		my @file_names = split(',', $ibams);		
		foreach my $bam_file (@file_names) {
			#$bam_file =~ /.+\/(.+?)\..+/;
			$bam_file =~ /([^\/]+)\..+/;
			my $prefix = $1;
			my $ifq = $odir."/$prefix.1.fq";
			my $ifq2 = $odir."/$prefix.2.fq";
			if ($ipfqs eq "") { $ipfqs = join (',', $ifq, $ifq2); }
			else { $ipfqs = join (',', $ipfqs, $ifq, $ifq2); }
			
			# convert the BAM file to fastq format 
			bam2fastq ($bam_file, $ifq, $ifq2);
		} 
		return $ipfqs;
	} # generate_input_fastq
	
	#############################################################################
	# convert each the BAM file to fastq format, the resulting file contains specified label $label 
	sub bam2fastq {
		my $bam_file = shift;
		my $ofq = shift;
		my $ofq2 = shift;
	
		# set NON_PF=True to use all reads instead of just PF reads
		
		print "\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=SILENT INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2\n\n";
		system ("\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=SILENT INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2");  	
	} # bam2fastq

	#############################################################################
	# given input list of fastq files generate output filenames for mvicuna
	sub generate_mvicuna_output_filenames {
		my $ref_drm_opfqs = shift;
		my $ref_prm_opfqs = shift;
		my $ref_prm_ofa = shift;
		my $ipfqs = shift;
		my $odir = shift;
		
		my @input_files = split (',', $ipfqs);
		foreach my $name (@input_files) {
			$name =~ /([^\/]+)\..+/;
			my $prefix = $1;
			if ($$ref_drm_opfqs eq "") { 
				$$ref_drm_opfqs = $odir."/$prefix.drm.fq";
				$$ref_prm_opfqs = "$odir/$prefix.drm.prm.1.fq,$odir/$prefix.drm.prm.2.fq";
				$$ref_prm_ofa .= "$odir/$prefix.drm.prm.fa"; 		
			} else { 
				$$ref_drm_opfqs .= ",$odir/$prefix.drm.fq"; 
			}			
		}
	} # generate_mvicuna_output_filenames
 