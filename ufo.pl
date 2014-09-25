# @brief
# Prerequisite: 
# 1) Perl 5.10 +
# 2) gcc-4.7.2 or above
# 3) set script path in ufo.pl


#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;
#use Digest::MD5 qw(md5);
 
	#--------------------------------------------------------------------------------------------------
	# set dependent programs/scripts/data path
	# system 
	my $perl_pl = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/perl_5.10.1/bin/perl";

	# dependent programs 
	my $sam2fastq_jar = "/seq/software/picard/current/bin/SamToFastq.jar";
	#my $trimmomatic_jar = "/seq/viral/analysis/xyang/external_programs/Trimmomatic-0.29/trimmomatic-0.29.jar";
	#my $trimref = "/seq/viral/analysis/xyang/FUO/DB/TruSeq3-PE.fa"; # Trimmomatic dependent	
	#my $bowtie_c = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/bowtie_0.12.9/bowtie";
	my $bwa_c = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/bwa_0.7.4/bwa";
	my $samtools_c = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/samtools/samtools_0.1.19/bin/samtools";
	my $runmosaik_pl = "/seq/viral/analysis/xyang/scripts/runMosaik.pl";
	my $dot_c = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/graphviz_2.26.3/bin/dot";
	
	# our package
	my $trimref = "/seq/viral/analysis/xyang/FUO/DB/adaptor_primer.fasta";
	#my $duprm_pl = "/seq/viral/analysis/xyang/FUO/scripts/duprm_by_sort.pl";
	#my $duprm_c = "/seq/viral/analysis/xyang/programs/DupRm/bin/duprm"; 
	my $mvicuna_c = "/seq/viral/analysis/xyang/programs/M-Vicuna/bin/mvicuna";
	my $cleanhits_pl = "/seq/viral/analysis/xyang/FUO/scripts/clean_hits.pl";
	my $rmhits_pl = "/seq/viral/analysis/xyang/FUO/scripts/rm_read_by_name.pl";
	my $rankabund_c = "/seq/viral/analysis/xyang/FUO/scripts//RankAbundance/bin/rankabund";
	#--------------------------------------------------------------------------------------------------

	my %option = (
		h           			=> '',
		quiet 					=> '',
		
		# I/O setting 
		ibams					=> '', # input comma separated BAM files, where each denotes a paired fq files
		ipfqs 					=> '', # input paired comma separated fastq files 
		isfqs					=> '', # input comma separated single fastq files
		odir					=> '', # output directory
		oprf					=> '', # output file prefix
		
		# performance setting 
		p 						=> 16,		 # performance control; num of cores to use
		batch					=> 1000000,   # performance control; batch of reads/read-pairs to be loaded in memory

		# preprocessing switch 
		no_preproc				=> '',

		# low complexity read setting 
		lc_n					=> '', # fraction of base n/N
		lc_mono					=> '', # fraction of mono-nt
		lc_di					=> '', # raction of any combination of 2 bases
		
		# duplicate removal setting
		drm_perc_sim			=> 95, # min perc similarity
		
		# trim setting 
		trm_vecfa 				=> '', # input vector adaptor primer fasta file 
		trm_min_match 			=> 9, # min match to trim 
		trm_min_rlen 			=> 70, # min read len to keep a read post-trimming
		trm_min_q				=> 2,  # min quality value 

		# de-contaminate using Mosaik, e.g. against ribosomal RNA
		dc_ref_fl				=> '', # reference sequence filelist used for de-contamination
		dc_mosaik_opt			=> "-m unique -act 20",
		
		# RankAbundance setting 
		rb_nodename				=> '',
		rb_tree					=> '',
		rb_skeleton				=> '',
		
		# aligner setting
		aln_only				=> '', # only to generate alignment result		
		aln_db_dir				=> '', # alignment database: every file in this folder will be aligned against	
		aln_db_fl				=> '', # file list to be aligned against
		aln_fraglen				=> 1000, # max fragment length
		aln_min_perc_span		=> 95, # min perc alignment region 
		aln_min_perc_sim 		=> 92,  # min perc similarity 		
	);

	GetOptions(
		"h"						=> \$option{h},
		"quiet"					=> \$option{quiet},
		
		# I/O setting
		"ibams=s"				=> \$option{ibams}, 	
		"ipfqs=s" 				=> \$option{ipfqs},	
		"isfqs=s"				=> \$option{isfqs},
		"odir=s"				=> \$option{odir}, 
		"oprf=s" 				=> \$option{oprf},
		
		# performance setting
		"p=i" 					=> \$option{p},  #number of processors to use
		"batch=i"				=> \$option{batch},   

		# switch for preprocessing, lc, duplrm, trim will not be carried out 
		"no_preproc"			=> \$option{no_preproc},
		
		# low complexity read setting
		"lc_n=f"				=> \$option{lc_n}, # fraction of base n/N
		"lc_mono=f"				=> \$option{lc_mono}, # fraction of mono-nt
		"lc_di=f"				=> \$option{lc_di}, # fraction of any combination of 2 bases 	
		
		# duplicate removal setting 
		"drm_perc_sim=i"		=> \$option{drm_perc_sim}, 
		
		# trim setting 	
		"trm_vecfa=s" 			=> \$option{trm_vecfa}, # input vector adaptor primer fasta file 
		"trm_min_match=i"  		=> \$option{trm_min_match}, # min match to trim 
		"trm_min_rlen=i" 		=> \$option{trm_min_rlen}, # min read len to keep a read post-trimming
		"trm_min_q=i" 			=> \$option{trm_min_q}, # min quality value 

		# de-contaminate using Mosaik, e.g. against ribosomal RNA
		"dc_ref_fl=s"			=> \$option{dc_ref_fl}, # reference sequences used for de-contamination
		"dc_mosaik_opt=s"		=> \$option{dc_mosaik_opt},
		
		# RankAbundance setting 
		"rb_skeleton=s" 			=> \$option{rb_skeleton}, # skeleton file
		"rb_nodename=s"			=> \$option{rb_nodename},
		"rb_tree=s"				=> \$option{rb_tree},
		
		# aligner setting
		"aln_only"				=> \$option{aln_only},
		"aln_db_dir=s"			=> \$option{aln_db_dir}, 
		"aln_db_fl=s" 			=> \$option{aln_db_fl},
		"aln_fraglen=i"			=> \$option{aln_fraglen},
		"aln_min_perc_span=i"	=> \$option{aln_min_perc_span}, # min perc alignment region 
		"aln_min_perc_sim=i" 	=> \$option{aln_min_perc_sim},  # min perc similarity 		

	) || printHelp (); 

	if ($option{h}) { printHelp();}

	# ----------------------------------------------- get input parameters -----------------------------------------------------------
	my $ipfqs = $option{ipfqs};
	my $isfqs = $option{isfqs};
	my $odir = $option{odir};
	my $oprf = $option{oprf};
	my $quiet_option = "";

	if ($option{quiet}) { $quiet_option = " -quiet ";  }

	unless ($option{ibams} || $ipfqs || $isfqs) { print "\n\tERR: one of -ibams, -ipfqs and -isfqs should be specified\n"; printHelp ();}
	unless ($odir) { print "\n\tERR:  -odir should be specified\n"; printHelp (); }
	unless ($oprf) { print "\n\tERR:  -oprf should be specified\n"; printHelp (); }
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	unless ($option{quiet}) { print "\n\nTime (min, hr, day) = ($min, $hour, $mday)\n\n"; }

	# --------------------- convert each bam file to paired fastq files, and record in $ipfqs --------------------------------
	
	if ($option{ibams}) { 
		if ($ipfqs) { $ipfqs .= ",".generate_input_fastq ($option{ibams}, $odir); }
		else { $ipfqs = generate_input_fastq ($option{ibams}, $odir); }
	} 
	
	# ----------------------------------------------------------------------------------------------------------------------------------------
	# - preprocessing via mvicuna: duplicate & low complexity read removal & trimming & paired read merging --
	# ----------------------------------------------------------------------------------------------------------------------------------------
	
	unless ($option{no_preproc}) {
		
		# ---------- prepare mvicuna output paired fastq files -----------
		my $drm_op;
		my $trm_op;
		my $trm_os;
		my $prm_op;
		my $prm_os;
		my $mvicuna_op = "$odir/$oprf.final.1.fq,$odir/$oprf.final.2.fq";
		my $mvicuna_os = "$odir/$oprf.final.s.fq";

		prep_mvicuna_ofilenames (\$drm_op, \$trm_op, \$trm_os, \$prm_op, \$prm_os, $odir."/$oprf");

		# ---------------- run mvicuna for preprocessing ------------------
		
		unless ($option{quiet}) { print "\n[[[ M-Vicuna: -tasks DupRm,Trim,PairedReadMerge]]]\n\n"; }
		
		my $lc_option;
		if ($option{lc_n}) {	$lc_option .= " -lc_n $option{lc_n} "; }
		if ($option{lc_mono}) { $lc_option .= " -lc_mono $option{lc_mono} "; }
		if ($option{lc_di}) { $lc_option .= " -lc_di $option{lc_di} "; }
	
		my $drm_option;
		if ($option{drm_perc_sim}) { $drm_option .= " -drm_perc_sim $option{drm_perc_sim} "; }
		my $trim_option;
		if ($option{trm_vecfa}) { $trim_option .= " -trm_vecfa $option{trm_vecfa} "; }
		if ($option{trm_min_match}) {	$trim_option .= " -trm_min_match $option{trm_min_match} "; }
		if ($option{trm_min_rlen}) { $trim_option .= " -trm_min_rlen $option{trm_min_rlen} "; }
		if ($option{trm_min_q}) { $trim_option .= " -trm_q $option{trm_min_q} "; }
	
		system ("$mvicuna_c -batch $option{batch} -pthreads $option{p} -ipfq $ipfqs -opfq $mvicuna_op -osfq $mvicuna_os -tasks \\
			 DupRm,Trim,PairedReadMerge $drm_option -drm_op $drm_op $lc_option -trm_op $trm_op -trm_os \\
			 $trm_os  $trim_option -prm_op $prm_op -prm_os $prm_os") == 0 or die "mvicuna failed"; 

		# output files: a paired fastq files stored as comma separated string in [$mvicuna_op] and a single-end fastq file [$mvicuna_os];

		$ipfqs = $mvicuna_op ;
		if ($isfqs) {	$isfqs .= ",".$mvicuna_os ;	}
		else { $isfqs = $mvicuna_os; }
		
	}  #	unless ($option{no_preproc}) 
 
  	# ----------------------------------------------------------------------------------------------------------------------------------------
	# de-contaminate using Mosaik, e.g. against ribosomal RNA, BWA aln or mem has very low sensitivity
	# ----------------------------------------------------------------------------------------------------------------------------------------
	unless ($option{quiet}) { print "\n[[[ De-contaminate rRNA, UniVec ]]]\n\n"; }
	
	# generate query file list 
	my @query_list;
	my $num_paired_end_fqs = scalar (split (',', $ipfqs));
	if ($isfqs) { 
		if ($num_paired_end_fqs > 0) { @query_list = split (',', $ipfqs.",$isfqs"); }
		else { @query_list = split (',', $isfqs); }
	} else { @query_list = split (',', $ipfqs); }
	unless ($option{quiet}) { print "\n\t".@query_list." query files IDed\n"; }
	
	my @dc_fl = split (',', $option{dc_ref_fl});
	my $sz_dc_fl = @dc_fl;
	my $sz_query_fl = @query_list;

	if ($sz_dc_fl > 0 && $sz_query_fl > 0) {
	
	 	my $dc_result_file = $odir."/$oprf.dc.txt";
		
		for (my $i = 0; $i < $sz_dc_fl; ++ $i) { # reference file
			for (my $j = 0; $j < $sz_query_fl; ++ $j) {	# query file	
				my $opath = $odir."/$oprf.dc.$i.$j.aln";
				system("perl $runmosaik_pl -p $option{p} -fq $query_list[$j] -ref $dc_fl[$i] -o $opath $option{dc_mosaik_opt}");
				my $rslt_bam = $opath.".sorted.bam";
				if (! -e $rslt_bam) {
					print "[ERR] Mosaik alignment failed to generate $rslt_bam\n";
					exit;
				} else {
					my $is_second_pair = 0;
					if (($j < $num_paired_end_fqs) && ($j % 2 == 1)) { $is_second_pair = 1; }
					
					# system ("$samtools_c view $rslt_bam | cut -f 1,3,4 | awk '{\$4=$is_second_pair; OFS=\"\t\"; print}' > $dc_result_file");	
	
					if (-e $dc_result_file) { system ("rm $dc_result_file"); }		
					system ("$samtools_c view $rslt_bam | cut -f 1 > $dc_result_file");	# only need read names

				}	
				`rm $rslt_bam`;	
				
				# remove hits from query file 
				my $tmp_query = $query_list[$j].".clean.fq";
				if (-e $tmp_query) { system ("rm $tmp_query"); } 
 				system ("$perl_pl $rmhits_pl -itsv $dc_result_file -ifq $query_list[$j] -ofq $tmp_query $quiet_option");
				system ("mv $tmp_query $query_list[$j]");
			}
		}		
	
		# fix read pairs
		if ($num_paired_end_fqs > 0) {
			unless ($option{quiet}) { print "\n\tFix Read Pairs in $ipfqs\n\n"; }
			fix_read_pairs ("$perl_pl $rmhits_pl", $ipfqs, "", $quiet_option);  # NOTE: if only one end has hits, whole read-pair is removed
		}
		
		# remove tmp files
		system ("rm $dc_result_file");
		
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
		unless ($option{quiet}) { print "\n\nTime (min, hr, day) = ($min, $hour, $mday)\n\n";}
	}

 	# ----------------------------------------------------------------------------------------------------------------------------------------
	# RankAbundance
	# ----------------------------------------------------------------------------------------------------------------------------------------
	unless ($option{quiet}) { print "\n[[[ RankAbundance ]]]\n\n"; }
	
	my $ipfqs_option = "";
	if ($ipfqs) { $ipfqs_option = "-ipfqs $ipfqs"; }
	system ("$rankabund_c -p $option{p} $ipfqs_option -isfqs $isfqs -iskeleton $option{rb_skeleton} -inodename $option{rb_nodename} \\
			-itree $option{rb_tree} -odir $option{odir} -oprefix $option{oprf}");
	
	#outputs taxa_tree_$option{oprf}.txt, taxa_rank_$option{oprf}.txt, taxa_rank_$option{oprf}.100.txt, spurious_taxa_$option{oprf}.txt
	
	unless ($option{quiet}) { print "\t$dot_c -Tpdf $odir/taxa_tree_$option{oprf}.txt -o $odir/taxa_tree_$option{oprf}.pdf\n"; }
	unless ($option{quiet}) { print "\t$dot_c -Tpdf $odir/taxa_tree_$option{oprf}.100.txt -o $odir/taxa_tree_$option{oprf}.100.pdf\n"; }

	system ("$dot_c -Tpdf $odir/taxa_tree_$option{oprf}.txt -o $odir/taxa_tree_$option{oprf}.pdf");
	system ("$dot_c -Tpdf $odir/taxa_tree_$option{oprf}.100.txt -o $odir/taxa_tree_$option{oprf}.100.pdf");
	system ("$dot_c -Tpng $odir/taxa_tree_$option{oprf}.txt -o $odir/taxa_tree_$option{oprf}.png");
	system ("$dot_c -Tpng $odir/taxa_tree_$option{oprf}.100.txt -o $odir/taxa_tree_$option{oprf}.100.png");
	
	print "debug done\n"; exit;
	# ----------------------------------------------------------------------------------------------------------------------------------------
	# ------------- All-read vs all database alignment via BWA 
	# ------------- Generate alignment results in format of [read_name]\t[ref_name]\t[hit_pos]\t[prj], where prj is the 
	# 		   projection using to distinguish two reads from a read pair
	# ------------- Generate Cleaned paired end files and singleton files
	# ----------------------------------------------------------------------------------------------------------------------------------------
	unless ($option{quiet}) { print "\n[[[ ALIGNMENT ]]]\n\n"; }
	
	# --- DB preparation: obtain all bwt index files from -aln_db_dir or/and -aln_db_fl ---
	
	my @db_list;
	if (-d $option{aln_db_dir}) { # reading from folder
		opendir (MYDIR, $option{aln_db_dir});
		my @dir_files = readdir (MYDIR);
		
		closedir(DH);
		
		my %bwa_index_files;
		foreach my $file (@dir_files) {
		
			next if($file =~ /^\.$/ || $file =~ /^\.\.$/ || $file =~ /pr\./); 	# skip . , .. and pr (protein)
			#$file =~ /[^\/][^.]([^.]*)/; # get file name a from path/a.b.c
			$file =~ /(.*\/)*(.*)\./;  # get $2=a.b from path/a.b.c

			unless (exists $bwa_index_files{$2}) {
				push (@db_list, $option{aln_db_dir}."/".$2);
				$bwa_index_files{$2} = 1;
			}
		}
	}

	if ($option{aln_db_fl}) { # reading from file list
		my $sz_db_list = @db_list;
		if ($sz_db_list > 0) { push (@db_list, split (',', $option{aln_db_fl}));	}
		else { @db_list = split (',', $option{aln_db_fl}); }		
	}
	
	# generate alignment tasks 
	my @aln_task_tuples; 
	my $idx_query = 0;
	foreach my $query (@query_list) {
		$query =~ /(.*\/)*(.*)\./;  # get $2=a.b from path/a.b.c
		my $query_prefix = $2;		
 		foreach my $db (@db_list) {
			$db =~ /(.*\/)*(.*)/;  # get $2=a.b.c from path/a.b.c
			my $db_prefix = $2;
			my $hitfile = $odir."/".$oprf."-$query_prefix"."-$db_prefix.txt";
			if (($idx_query < $num_paired_end_fqs) && ($idx_query % 2 == 1)) { # second pair
				push (@aln_task_tuples,  ($query, $db, $hitfile, 1));   
			} else { push (@aln_task_tuples,  ($query, $db, $hitfile, 0));  } 
		}
		++ $idx_query;
	}
 

	# ------------- carry out alignment for individual input fastq file and merge all alignment results --------------------
	
	unless ($option{quiet}) { print "\n\t".(@aln_task_tuples/4)." BWA alignments to be carried out\n\n"; }

	my $num_tasks = scalar (@aln_task_tuples) / 4;

	if ($num_tasks > 0) {
		my $aln_result = $odir."/$oprf.bwa_aln_rslt.txt";
		if (-e $aln_result) { system ("rm $aln_result"); }
		my $min_score = $option{trm_min_rlen} * $option{aln_min_perc_span} * $option{aln_min_perc_sim} / ( 100 * 100 );	
		my $fm = Parallel::ForkManager->new ($option{p});
		for (my $i = 0; $i < $num_tasks; ++ $i) {  

			my $pid = $fm->start and next;
	
			my $bwa_cmd_all = "$bwa_c mem -v 1 -a -k 27 -w 6 -T $min_score $aln_task_tuples[4*$i + 1] $aln_task_tuples[4* $i]"; # -w band, -T min_score					
			my $bwa_cmd = "$bwa_c mem -v 1 -k 27 -r 2 -w 6 -T $min_score $aln_task_tuples[4*$i + 1] $aln_task_tuples[4* $i]"; # -w band, -T min_score					
			$bwa_cmd = "$bwa_c aln $aln_task_tuples[4*$i + 1] $aln_task_tuples[4* $i]"; # -w band, -T min_score					

			#"$bwa_c mem -a -t $option{p} -k 23 -w 10 -T $min_score $aln_task_tuples[4*$i + 1] $aln_task_tuples[4* $i]"; # -w band, -T min_score			
		
			unless ($option{quiet}) { print "\n\t$i : align $aln_task_tuples[4*$i] against $aln_task_tuples[4*$i + 1]\n"; }	
		
			run_bwa_mem ($bwa_cmd, $option{aln_min_perc_span}, $option{aln_min_perc_sim}, $aln_task_tuples[4*$i + 2], $aln_task_tuples[4*$i + 3]) ;		
		
			$fm->finish;
		}
		$fm->wait_all_children;	

		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
		unless ($option{quiet}) { print "\n\nTime (min, hr, day) = ($min, $hour, $mday)\n\n"; }

		# --------------------------- merge results ------------------------------------
		for (my $i = 0; $i < $num_tasks; ++ $i) {  	
			system ("cat $aln_task_tuples[4*$i + 2] >> $aln_result");  # merge alignment results 
			system ("rm $aln_task_tuples[4*$i + 2]");		
		}
	
		# -------------------- sort alignment result by (read name, refname) -----------------------
		unless ($option{quiet}) { print "\n\t Sort alignment result\n"; }
		my $all_aln_tmp = $odir."/$oprf.bwa_aln_rslt.sort.txt";
		system ("sort -k 1,2 $aln_result > $all_aln_tmp");
		system ("mv $all_aln_tmp $aln_result");
	
		if ($option{aln_only}) { 
			print "\tcounting alignments in $aln_result\n";		
			count_aln ($aln_result, $option{aln_fraglen}); # counting aligned reads, and fragments
			print "Done\n"; exit; 
		}
	
		# -------------------- analyze aln results and clean hits ------------------------------
		unless ($option{quiet}) { 
			print "\n\t[ Clean alignment result ]\n"; 
			print "\n\t\t$perl_pl $cleanhits_pl -maxflen $option{aln_fraglen} -i $aln_result -o $all_aln_tmp -p $option{p}\n";
		}	
		system ("$perl_pl $cleanhits_pl -maxflen $option{aln_fraglen} -i $aln_result -o $all_aln_tmp -p $option{p}");
		system ("mv $all_aln_tmp $aln_result");
	
		# ---------------- clean out reads with hits from input files --------------------------
		unless ($option{quiet}) { print "\n\t[ Clean hits from input files ]\n"; }
		foreach my $query (@query_list) { 
			my $tmpfile = $query.".clean.fq";
			if (-e $tmpfile) { system ("rm $tmpfile"); } 
			unless ($option{quiet}) { print "\n\tperl $rmhits_pl -itsv $aln_result -ifq $query -ofq $tmpfile $quiet_option\n"; }		
			system ("$perl_pl $rmhits_pl -itsv $aln_result -ifq $query -ofq $tmpfile $quiet_option");
			system ("mv $tmpfile $query");		
		}
		# -------------- fix read pairs -----------------------------
		if ($num_paired_end_fqs > 0) {
			unless ($option{quiet}) { print "\n\t[[[ Fix Read Pairs in $ipfqs ]]]\n\n"; }
			fix_read_pairs ("$perl_pl $rmhits_pl", $ipfqs, "", $quiet_option);  # NOTE: if only one end has hits, whole read-pair is removed
		}
		
		($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
		unless ($option{quiet}) { print "\n\nTime (min, hr, day) = ($min, $hour, $mday)\n\n";}
	
	} # if ($num_tasks > 0)

	print "\ndebug done\n";
	exit;


	
	# --- generate { readname -> read_cnt } hash table from alignment result [todo] the count needs to be associate with phylogenetic tree ---	
	#my %rname_cnt;
	#gen_readcnt_hash (\%rname_cnt, $aln_result);
	
		
	# --- loop: get hits to current database then remove these hits from further consideration ----------
#	foreach my $db (@db_list) {	
#		$db =~ /(.*\/)*(.*)\./;  # get $2=a.b from path/a.b.c
#		my $db_prefix = $2;
#		my $hitfile = $odir."/".$oprf.".$db_prefix.hit.txt";
#		if (-e $hitfile) { system ("rm $hitfile"); }
#		my $input_fqs = $mvicuna_op.",$mvicuna_os";
#		my @query_fqs = split (',', $input_fqs);
#		foreach my $query_fq (@query_fqs) {
#			my $min_score = $option{trm_min_rlen} * 0.7;	
#			my $bwa_cmd = "$bwa_c mem -t $option{p} -k 23 -w 10 -T $min_score $db $query_fq"; # -w band, -T min_score			
#			run_bwa_mem ($bwa_cmd, $option{aln_min_perc_span}, $option{aln_min_perc_sim}, $hitfile);		
#			my $tmp_ofq = $query_fq.".clean.txt";			
#			if (-e $tmp_ofq) { system ("rm $tmp_ofq"); } 
#			unless ($option{quiet}) { print "\n\tperl $rmhits_pl -itsv $hitfile -ifq $query_fq -ofq $tmp_ofq $quiet_option\n"; }		
#			system ("$perl_pl $rmhits_pl -itsv $hitfile -ifq $query_fq -ofq $tmp_ofq $quiet_option");
#			system ("mv $tmp_ofq $query_fq");
#		}		
#	} # foreach my $db

	# fix read pairs 
#	unless ($option{quiet}) { print "\n[[[ fix_read_pairs in $mvicuna_op ]]]\n\n"; }
#	fix_read_pairs ("$perl_pl $rmhits_pl", $mvicuna_op, "", $quiet_option);


	#################################################################################################################################
	# sub-functions 
	#############################################################################
 
	sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "usage: perl $0 {-ibams [bamfilelist] or -ipfqs [pairedFQfilelist]} -odir [odir/] -oprf [oprefix] \n\n";
	
		print "-h: help; for cmd options\n";
		print "-quiet: no screen output for programs called\n";
	
		print "\n*** I/O setting ***\n";		
		print "-ibams: comma separated BAM files; \n";
		print "-ipfqs: comma separated paired fq files\n";
		print "-isfqs: comma separated singleton fq files\n";
		print "-odir: output directory\n";
		print "-oprf: output file prefix\n";

		print "\n*** Performance setting ***\n";
		print "-p: default 16; number of threads\n";
		print "-batch: default 1 mil; number of reads/read-pairs to be in memory\n";

		print "\n*** Preprocessing: duplicate, low complexity removal and trim setting, used by mvicuna ***\n";
		print "-no_preproc: default unspecified, if specified, initial duplicatie and low complexity\n\t\t removal and trimming will not be done\n";
		print "-lc_n: default 40; percentage of base n/N\n";
		print "-lc_mono: default 50; percentage of mono-nt\n";
		print "-lc_di: default 80; percentage of any combination of 2 bases\n";						
		print "-drm_perc_sim: default 95; min percent similarity \n";
		print "-trm_vecfa: optional; input vector adaptor primer fasta file\n";
		print "-trm_min_match: default 9; min length of match required to trim \n";
		print "-trm_min_rlen: defaut 70; min read len to keep a read post-trimming\n";
		print "-trm_min_q: defaut 2 (phred score)\n";
	
		print "\n*** De-contamination using Mosaik, e.g. against ribosomal RNA ***\n";
		print "-dc_ref_fl: comma separated ref sequence file list\n";
		print "-dc_mosaik_opt: default '-m unique -act 20'; mosaik options\n";
		
		print "\n*** RankAbundance setting ***\n";
		print "-rb_skeleton: skeleton file\n";
		print "-rb_nodename: file of node names, e.g. names.clean.dmp\n";
		print "-rb_tree: file specifying taxonomy tree, e.g. nodes.dmp\n";
				
		print "\n*** Aligner setting ***\n";
		print "-aln_only: only perform and output alignment conforming specified criteria\n";
		print "-aln_db_dir: folder containing all BWA index files to be used\n";
		print "-aln_db_fl: comma separated BWA index files to be used\n";
		print "-aln_fraglen: default 1000, max fragment length\n";
		print "-aln_min_perc_span: default 95; min perc of alignment region\n";
		print "-aln_min_perc_sim: default 92; min perc similarity of alignment region\n"; 		
								
		print "\n----------------------------------------------------------------------------\n\n";
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
		
		print "\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=quiet INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2\n\n";
		system ("\n\tjava -jar $sam2fastq_jar QUIET=true MAX_RECORDS_IN_RAM=$option{batch} VALIDATION_STRINGENCY=quiet INPUT=$bam_file FASTQ=$ofq SECOND_END_FASTQ=$ofq2");  	
	} # bam2fastq

	#############################################################################
	# given input list of fastq files generate output filenames for mvicuna
	sub prep_mvicuna_ofilenames {
		my $ref_drm_opfqs = shift;
		my $ref_trm_opfqs = shift;
		my $ref_trm_os = shift;
		my $ref_prm_opfqs = shift;
		my $ref_prm_os = shift;
		my $opath = shift;

		$$ref_drm_opfqs = "$opath.drm.1.fq,$opath.drm.2.fq";
		$$ref_trm_opfqs = "$opath.trm.1.fq,$opath.trm.2.fq";
		$$ref_trm_os = "$opath.trm.s.fq";
		$$ref_prm_opfqs = "$opath.prm.1.fq,$opath.prm.2.fq";
		$$ref_prm_os = "$opath.prm.s.fq"; 		
	} # prep_mvicuna_ofilenames


	#############################################################################
	# run bowtie for a list of input fastq files 
	sub run_bowtie {
 		my $cmd = shift;
 		my $db = shift;
 		my $opt = shift;
 		my $inputs = shift;
 		my $output = shift;
	
		my @input_files = split (',', $inputs);
		foreach my $input_file (@input_files) {
			system ("$cmd $db $opt -q $input_file >> $output");
		}
  	} # run_bowtie

 	#############################################################################
 	# return the number of specific cigar character in a cigar string 
 	# sub-function called by run_bwa_mem 
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
	
	# run bwa mem for a list of input fastq files 
	# input 1) bwa_cmd, 2) min_perc_span, 3) min_perc_sim, 4) output file, 5) which end
	sub run_bwa_mem {
		my $bwa_cmd = shift;
		my $min_perc_span = shift;
		my $min_perc_sim = shift;
		my $ofile = shift;
		my $which_end = shift;
		
		open(OUTPUT, ">$ofile") or die $!;
		open(CMD, "$bwa_cmd |") or die $!;
		my @buffer;
		while (my $line = <CMD>){
			my @entries = split ('\t', $line);
			unless ($entries[0] =~ /^@/){
				if ($#entries + 1 < 11) {
					print "[error] line: $line has < 11 entries\n";
					exit;
				}
				my $query_name = $entries[0];			
				#my $flag = $entries[1];     # used to indicate if this is the first end or the last end 0x40: first  	
				my $hit_name = $entries[2];				
				my $hit_pos = $entries[3];
				my $cigar = $entries[5];
				my $query_seq = $entries[9];	
				my $num_match = get_cigar_char_cnt ($cigar, 'M');
				next if $num_match == 0;
				my $query_len = length ($query_seq);	
				my $query_seq_len = $query_len + get_cigar_char_cnt ($cigar, 'H');	# raw sequence input length
				my $span_len = $query_len - get_cigar_char_cnt ($cigar, 'S'); # minus soft-clip
				my $perc_span = $span_len * 100 / $query_seq_len;
				my $perc_sim =  $num_match * 100 / get_cigar_char_cnt ($cigar, "[M,I,D,N,=,X]");
				if (($perc_span >= $min_perc_span) && ($perc_sim >= $min_perc_sim)) {

					#my $str = substr( md5($query_seq), 0, 4 );	
					#my $prj = unpack('L', $str) % 1000000;
														
					if ($hit_name =~ /gi\|(\d+)\|/) {
						my $gi = $1;				
					#	print OUTPUT $query_name."\t".$gi."\t"."$perc_span\t$perc_sim\n";		
					#	print OUTPUT "$query_name\t\t$gi\t\t$hit_pos\t\t$prj\t\t$perc_span\t\t$perc_sim\n";
					#	print OUTPUT "$query_name\t$gi\t$hit_pos\t$which_end\n";
						push (@buffer, "$query_name\t$gi\t$hit_pos\t$which_end\n");
					} else { # no gi number 
					#	print OUTPUT "$query_name\t$hit_name\t$hit_pos\t$which_end\n";
						push (@buffer, "$query_name\t$hit_name\t$hit_pos\t$which_end\n");						
					}
					# output
					if (scalar (@buffer) > 100000) {
						foreach my $elem (@buffer) {
							print OUTPUT $elem;							 
						}
						@buffer=();
					}
				}		
			} # unless 
	
		} # while 
		foreach my $elem (@buffer) {
			print OUTPUT $elem;
		}
		@buffer=();
		close (CMD);
		close (OUTPUT);
	} # run_bwa_mem

 	#############################################################################
 	# Given bwa alignment record file: { [rname] [ref_gi] [aln_span] [aln_sim] } sorted by ( [rname], [ref_gi] )
 	# generate hash table recording  { rname -> count }	
	sub gen_readcnt_hash {
		my $ref_hash = shift;
		my $ifile = shift;
		open (TSV, $ifile) or die "[ERR] gen_hash_readcnt, cannot open $ifile\n";
		my $prev_rname;
		my $prev_ref;
		my $cnt = 1; 
		while (my $line = <TSV>) {
			my @entries = split ('\t', $line);
			if ($prev_rname) {
				if ($entries[0] eq $prev_rname) { ++ $cnt unless ($entries[0] eq $prev_ref); } 				
				else {
					$ref_hash->{$prev_rname} = $cnt;
					$cnt = 1;
					$prev_rname = $entries[0];
					$prev_ref = $entries[1];
				}
			} else { # first entry 
				$prev_rname = $entries[0];
				$prev_ref = $entries[1];				
			}
		}
		$ref_hash->{$prev_rname} = $cnt;		
		close (TSV);
	} #  gen_readcnt_has

 	#############################################################################			
 	# Given bwa alignment record file {[rname] [refname] [startpos] [is_first_end]}
 	#		sorted by rname and refname and max fragment length
 	# calculate # aligned reads, and # aligned fragments
	sub	count_aln {
		my $ifile = shift;
		my $fraglen = shift;
		open (TSV, $ifile) or die "[ERR] count_aln, cannot open $ifile\n";
		my @prev_entry;
		my $is_1_aligned = 0;
		my $is_2_aligned = 0;
		my $is_both_aligned = 0;
		my $single_cnt = 0;
		my $pair_cnt = 0;
		my $debug_cnt = 0;
		while (my $line = <TSV>) {
			
			my @entry = split ('\t', $line);
			chomp (@entry);
  	
#			++ $debug_cnt;
#			if ($debug_cnt == 12) {
#				print $single_cnt.", ".$pair_cnt."\n";
#				exit;
#			}

			if (scalar (@prev_entry) > 0) {
				if ($entry[0] eq $prev_entry[0]) { 
					if ($is_both_aligned != 1) {
					
						#if ($prev_entry[1] == 224589813) {
						#	print "here\n";
						#	print Dumper(@prev_entry, @entry);			
						#}
						
						if ( ($entry[1] eq $prev_entry[1]) && # same ref
							 (abs($entry[2] - $prev_entry[2]) <= $fraglen) && # fragment constraint
							 ($entry[3] != $prev_entry[3]) )  # different ends
						{ $is_both_aligned = 1; }				
						else {
							if ($entry[3] == 0) { $is_1_aligned = 1; }
							else { $is_2_aligned = 1; }
						}
					} 
				} else {
					if ($is_both_aligned == 1) { ++ $pair_cnt; }
					else {	
						if ($is_1_aligned == 1) { ++ $single_cnt; }
						if ($is_2_aligned == 1) { ++ $single_cnt; }
					} 
					$is_1_aligned = 0;
					$is_2_aligned = 0;
					$is_both_aligned = 0;					
					if ($entry[3] == 0) { $is_1_aligned = 1;}
					else { $is_2_aligned = 1; }			
					
#					print Dumper(@prev_entry, @entry);		
#					print $single_cnt.", ".$pair_cnt."\n";
#					++ $debug_cnt;
#					if ($debug_cnt == 20) {
#						exit;
#					}
				}
				@prev_entry = @entry;				
			} else { # first entry 
				@prev_entry = @entry;
				if ($prev_entry[3] == 0) { $is_1_aligned = 1;}
				else { $is_2_aligned = 1; }
			}
		}				
		close (TSV);
		if ($is_both_aligned == 1) { ++ $pair_cnt; } 
		else {	
			if ($is_1_aligned == 1) { ++ $single_cnt; }
			if ($is_2_aligned == 1) { ++ $single_cnt; }
		} 
	
		print "\t\tTotal aligned singletons, pairs, and total: ";
		print $single_cnt.", ".$pair_cnt.", ".($single_cnt + 2*$pair_cnt)."\n";
		 
	}
 	#############################################################################
 	# Given a paired fastq input files, remove singletons from both files;
 	# ADD these singletons to a single fq file $sfq if it is not empty
 	sub fix_read_pairs {
 		my $rmhits_pl = shift;
	 	my $pfq = shift;
	 	my $sfq = shift;
	 	my $quiet_option = shift;
	 	
	 	# obtain two tsv files storing singleton read names in fastq file pair
		my @tmp_tsv_files = get_singleton_reads ($pfq, $quiet_option); 		
		my @ipfq_files = split (',', $pfq);
		my @tmp_opfq_files;
		foreach my $ifq (@ipfq_files) { push (@tmp_opfq_files, $ifq.".clean.fq"); }
		
		for (my $i = 0; $i < 2; ++ $i) {
			if ($sfq) {	# add singletons to $sfq
				print "\t$rmhits_pl -itsv $tmp_tsv_files[$i] -ifq $ipfq_files[$i] -ofq $tmp_opfq_files[$i] -ofqrm $sfq $quiet_option\n";
				system ("$rmhits_pl -itsv $tmp_tsv_files[$i] -ifq $ipfq_files[$i] -ofq $tmp_opfq_files[$i] -ofqrm $sfq $quiet_option");	 	
			} else { # discard singletons
				print "\t$rmhits_pl -itsv $tmp_tsv_files[$i] -ifq $ipfq_files[$i] -ofq $tmp_opfq_files[$i] $quiet_option\n";
				system ("$rmhits_pl -itsv $tmp_tsv_files[$i] -ifq $ipfq_files[$i] -ofq $tmp_opfq_files[$i] $quiet_option");	 	
			}
			system ("rm $tmp_tsv_files[$i]");
			system ("mv $tmp_opfq_files[$i] $ipfq_files[$i]");
		}
		
 	} # fix_read_pairs
 	
 	# Given a comma separated input fastq files, identify singletons in both files
 	# and output in 2 tsv files where readname is the first column
 	# return comma concatenated 2 tsv filenames
 	sub	get_singleton_reads {
 		my $ipfq = shift;
 		my $quiet_option = shift;
 		my @otsv;
 		 		
 		my @pairs = split (',', $ipfq);
 		if ($#pairs + 1 != 2) {
 			print "[error] sub func: get_singleton_reads() input paired fqs $ipfq size !=2\n";
 			exit;
 		}
 		
 		my %rnames;
 		my %rnames2; 		
		get_fq_read_names (\%rnames, $pairs[0]);
		get_fq_read_names (\%rnames2, $pairs[1]);		

		my $tsv_output = $pairs[0].".tsv";
		push (@otsv, $tsv_output);
 		open my $fh_otsv, ">$tsv_output" or die $!; # output tsv
 		$tsv_output = $pairs[1].".tsv";
		push (@otsv, $tsv_output); 		
		open my $fh2_otsv, ">$tsv_output" or die $!; # output tsv
	
		my $num_s = 0;
		my $num_s2 = 0;			
		foreach my $key (keys %rnames) {
			unless (exists $rnames2{$key}) { # not found 
				print $fh_otsv $key."\n";
				++ $num_s;
			}
		}		
		foreach my $key (keys %rnames2) {
			unless (exists $rnames{$key}) { # not found 
				print $fh2_otsv $key."\n";
				++ $num_s2;
			}
		}		
		close ($fh_otsv);
		close ($fh2_otsv);
		
		my $num_keys = scalar keys %rnames;
		unless ($quiet_option) { print "\t$pairs[0]: num_reads, num_singletons = $num_keys, $num_s\n\n"; }
		$num_keys = scalar keys %rnames2;
		unless ($quiet_option) { print "\t$pairs[1]: num_reads, num_singletons = $num_keys, $num_s2\n\n"; }

		return @otsv;
		
	} # get_singleton_reads
	
	sub get_fq_read_names {
		my $ref_rname = shift;
		my $filename = shift;
		# go through fastq file
		open (FQ, $filename) or die "Cannot open $filename\n";
		while (my $line = <FQ>) {
			if ($line =~ /^\@/) {
				my $header = $line;
				chomp ($header);
				# rm first char '@'
				$header = substr($header, 1, length($header) - 1); 
				if ($header =~ /\/1$/ || $header =~ /\/2$/) {
					# rm /1 or /2
					$header = substr ($header, 0, length ($header) - 2); 
				}
				${$ref_rname} {$header} = 1;
			
				# skip the remaining lines 
				for (my $i = 0; $i < 3; ++ $i) { $line = <FQ>;}
			} else { die "\npotential error in the format of file $filename\n\n"; }
		} # while 		
	}

