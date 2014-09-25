
1. Denovo duplicate & low complexity read removal -- done about 23 mins on 8 cores. 20% dupl removal

2. Trimmomatic -- done till this step about half an hr. 

---------- Existing setting ----------

TrimmomaticPE: Started with arguments: -threads 8 -phred33 /idi/sabeti-data/stremlau/fastq/49CMiSeq.reads1.fastq /idi/sabeti-data/stremlau/fastq/49CMiSeq.reads2.fastq 49cmiseq.trim.p1.fq 49cmiseq.trim.s1.fq 49cmiseq.trim.p2.fq 49cmiseq.trim.s2.fq LEADING:25 TRAILING:25 SLIDINGWINDOW:4:25 MINLEN:70 ILLUMINACLIP:/seq/viral/analysis/xyang/FUO/DB/TruSeq3-PE.fa:2:30:10
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 17617141 Both Surviving: 7609520 (43.19%) Forward Only Surviving: 5219641 (29.63%) Reverse Only Surviving: 603552 (3.43%) Dropped: 4184428 (23.75%)

+++++++++++ New setting ++++++++++
TrimmomaticPE: Started with arguments: -threads 8 -phred33 /idi/sabeti-data/stremlau/fastq/49CMiSeq.reads1.fastq /idi/sabeti-data/stremlau/fastq/49CMiSeq.reads2.fastq 49cmiseq.trim.p1.fq 49cmiseq.trim.s1.fq 49cmiseq.trim.p2.fq 49cmiseq.trim.s2.fq LEADING:15 TRAILING:15 SLIDINGWINDOW:10:15 MINLEN:50 ILLUMINACLIP:/seq/viral/analysis/xyang/FUO/DB/TruSeq3-PE.fa:2:30:10
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 17617141 Both Surviving: 15007484 (85.19%) Forward Only Surviving: 2412439 (13.69%) Reverse Only Surviving: 116738 (0.66%) Dropped: 80480 (0.46%)
TrimmomaticPE: Completed successfully
599.470u 53.231s 1:48.06 604.0%	0+0k 0+0io 1pf+0w

3. Hit Analysis -- searching against known database ; till this step about 1 hr

a) bowtie alignment 
	bowtie DB/bowtie_idx_metagenomics_contaminants_v2 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.s.fq > tmp2/49cmiseq.hit.meta.txt
	bowtie DB/bowtie_idx_metagenomics_contaminants_v2 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p1.fq >> tmp2/49cmiseq.hit.meta.txt
	bowtie DB/bowtie_idx_metagenomics_contaminants_v2 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p2.fq >> tmp2/49cmiseq.hit.meta.txt
	bowtie DB/bowtie_idx_hg19 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.s.fq > tmp2/49cmiseq.hit.hg.txt
	bowtie DB/bowtie_idx_hg19 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p1.fq >> tmp2/49cmiseq.hit.hg.txt
	bowtie DB/bowtie_idx_hg19 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p2.fq >> tmp2/49cmiseq.hit.hg.txt
	bowtie DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.s.fq > tmp2/49cmiseq.hit.hrna.txt
	bowtie DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p1.fq >> tmp2/49cmiseq.hit.hrna.txt
	bowtie DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p2.fq >> tmp2/49cmiseq.hit.hrna.txt
	bowtie DB/bowtie_idx_nt_viruses -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.s.fq > tmp2/49cmiseq.hit.viral.txt
	bowtie DB/bowtie_idx_nt_viruses -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p1.fq >> tmp2/49cmiseq.hit.viral.txt
	bowtie DB/bowtie_idx_nt_viruses -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/49cmiseq.trim.lcrm.p2.fq >> tmp2/49cmiseq.hit.viral.txt


	bowtie DB/bowtie_idx_metagenomics_contaminants_v2 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo1.fq > tmp2/rabdo.hit.meta.txt
	bowtie DB/bowtie_idx_metagenomics_contaminants_v2 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo2.fq  >> tmp2/rabdo.hit.meta.txt
	bowtie DB/bowtie_idx_hg19 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo1.fq > tmp2/rabdo.hit.hg.txt
	bowtie DB/bowtie_idx_hg19 -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo2.fq >> tmp2/rabdo.hit.hg.txt
	bowtie DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo1.fq > tmp2/rabdo.hit.hrna.txt
	bowtie DB/bowtie_idx_ncRNA_transcripts_rRNA_mitRNA -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo2.fq >> tmp2/rabdo.hit.hrna.txt
	bowtie DB/bowtie_idx_nt_viruses -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo1.fq > tmp2/rabdo.hit.viral.txt
	bowtie DB/bowtie_idx_nt_viruses -p 12 -k 1 -v 2 --suppress 5,6,7,8 -q tmp2/rabdo2.fq >> tmp2/rabdo.hit.viral.txt

b) pairwise comparing the hit files, obtain shared elements (.shared.txt), summarize hits to each reference (.refcnts.txt) -- [cmp_hit_files.pl]
	databases: 1) 49cmiseq.hit.hg.txt 2) 49cmiseq.hit.hrna.txt 3) 49cmiseq.hit.meta.txt 4) 49cmiseq.hit.viral.txt

	-- output pairwisely number of common reads sharing between database hits
	49cmiseq.hit.hg.txt	49cmiseq.hit.hrna.txt	3134788
	49cmiseq.hit.hg.txt	49cmiseq.hit.meta.txt	2996380
	49cmiseq.hit.hg.txt	49cmiseq.hit.viral.txt	3532
	49cmiseq.hit.hrna.txt	49cmiseq.hit.meta.txt	3652761
	49cmiseq.hit.hrna.txt	49cmiseq.hit.viral.txt	3907
	49cmiseq.hit.meta.txt	49cmiseq.hit.viral.txt	352122

	-- generate file: .refcnts. listing for each ref in the db the num of read hits
	gi|392790|gb|U00220.1|U00220	100218	Human immunodeficiency virus type 1 Nef protein and neomycin phosphotransferase genes, complete cds
	gi|296556485|gb|M19921.2|HIVNL43	41912	Human immunodeficiency virus type 1, NY5/BRU (LAV-1) recombinant clone pNL4-3
	...
	
c) given read names ((tsv file: read name as the 1st field)), obtain actual reads, generate read statistics by parsing through a given set of fq files -- [get_hit_stats.pl] 
 