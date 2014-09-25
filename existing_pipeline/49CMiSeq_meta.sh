#-------- MEGAN ANALYSIS @ BROAD --------#
# identify filtering reads
# NO TRAILING / IN URL
use BLAST+
use Perl-5.10
use Samtools
use Python-2.7
use Java-1.6
use Bowtie
use BWA

# MAKE REQUIRED SUB-DIRECTORIES - DON'T DO THIS IF YOU ALREADY CREATED THESE WITH ANOTHER PIPELINE
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
bsub -o ~/log.txt -P sabeti_meta -G sabeti_labfolk "mkdir $url/_logs $url/_temp $url/_reads $url/_bams $url/_reports $url/_pileup $url/_meta"
done

#-------- BAM FILES AS INPUT --------#
# CONVERT BAM FILES TO FASTQ
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for suffix in bam
do
for in_url in /idi/sabeti-data/stremlau/a121026_LIB5/_bams
do
bsub -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.fq "java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$in_url/$reads.$suffix FASTQ=$url/_reads/$reads.reads1.fastq SECOND_END_FASTQ=$url/_reads/$reads.reads2.fastq VALIDATION_STRINGENCY=SILENT"
done
done
done
done

#-------- FASTQ & BAM FILES AS INPUT --------#
# GENERATE QUALITY REPORT
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for suffix in reads1.fastq reads2.fastq
do
bsub -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.qc "/idi/sabeti-scratch/kandersen/bin/fastqc/fastqc -f fastq $url/_reads/$reads.$suffix -o $url/_reports/"
done
done
done

# TRIM THE READS WITH TRIMMOMATIC
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
for suffix in fastq
do
for in_url in $url/_reads
do
for phred in -phred33 # Alternatively -phred64
do
for minimum_length in 70
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.tr "java -Xmx2g -classpath /idi/sabeti-scratch/kandersen/bin/trimmomatic-0.15/trim.jar org.usadellab.trimmomatic.TrimmomaticPE $phred $in_url/$reads.reads1.$suffix $in_url/$reads.reads2.$suffix $temp/$reads.trimmed.1.fastq $temp/$reads.reads1.trimmed_unpaired.fastq $temp/$reads.trimmed.2.fastq $temp/$reads.reads2.trimmed_unpaired.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:$minimum_length ILLUMINACLIP:/idi/sabeti-scratch/kandersen/references/contaminants/contaminants.fasta:2:40:12"
done
done
done
done
done
done
done
done
done
done

# BMTAGGER REMOVAL OF HUMAN TRANSCRIPTS, rRNA AND mitRNA
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in /idi/sabeti-scratch/kandersen/references/bmtagger/GRCh37.68_ncRNA-GRCh37.68_transcripts-HS_rRNA_mitRNA
do
bsub -R "rusage[mem=8]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bt -E "ls /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh" "/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b $db.bitmask -x $db.srprism -T $temp -q1 -1 $temp/$reads.trimmed.1.fastq -2 $temp/$reads.trimmed.2.fastq -o $temp/$reads.bmtagger.mrna"
done
done
done
done

# BMTAGGER REMOVAL OF HUMAN GENOMIC READS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in /idi/sabeti-scratch/kandersen/references/bmtagger/hg19
do
bsub -R "rusage[mem=8]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bt -E "ls /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh" "/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b $db.bitmask -x $db.srprism -T $temp -q1 -1 $temp/$reads.bmtagger.mrna.1.fastq -2 $temp/$reads.bmtagger.mrna.2.fastq -o $temp/$reads.bmtagger.hg19"
done
done
done
done

# BMTAGGER REMOVAL OF CONTAMINATING READS
# This will remove Lassa - if you want to keep, use /idi/sabeti-scratch/kandersen/references/bmtagger/metagenomics_contaminants_v2
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in /idi/sabeti-scratch/kandersen/references/bmtagger/metagenomics_contaminants_v2
do
bsub -R "rusage[mem=8]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bt -E "ls /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh" "/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b $db.bitmask -x $db.srprism -T $temp -q1 -1 $temp/$reads.bmtagger.hg19.1.fastq -2 $temp/$reads.bmtagger.hg19.2.fastq -o $temp/$reads.bmtagger.contaminants"
done
done
done
done

# REMOVE DUPLICATES
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p1 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -ns_max_n 1 -derep 1 -fastq $temp/$reads.bmtagger.contaminants.1.fastq -out_bad null -line_width 0 -out_good $temp/$reads.prinseq.1"
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p2 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -ns_max_n 1 -derep 1 -fastq $temp/$reads.bmtagger.contaminants.2.fastq -out_bad null -line_width 0 -out_good $temp/$reads.prinseq.2"
done
done
done
done
done
done

# FIX MATE-PAIR INFORMATION--NOTE NEW 'REGULAR EXPRESSION'
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.fm "/idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastqSeqs.pl -t -r '^@(\S+)\s[1|2]\S+$' -f1 $temp/$reads.prinseq.1.fastq -f2 $temp/$reads.prinseq.2.fastq -o $url/_reads/$reads.clean"
done
done
done
done



#-------- CONTIG-BASED ANALYSIS --------#
# DE NOVO ASSEMBLY USING METAVELVET - 1
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.v1 "/idi/sabeti-scratch/kandersen/bin/velvet/velveth $temp/$reads.velvet 23 -shortPaired -fastq -separate $url/_reads/$reads.clean.1.fastq $url/_reads/$reads.clean.2.fastq"
done
done
done
done
done
done

# DE NOVO ASSEMBLY USING METAVELVET - 2
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.v2 "/idi/sabeti-scratch/kandersen/bin/velvet/velvetg $temp/$reads.velvet -read_trkg yes -ins_length 380 -cov_cutoff auto"
done
done
done
done
done
done

# DE NOVO ASSEMBLY USING METAVELVET - 3
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
for contig_length in 200
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.v3 "/idi/sabeti-scratch/kandersen/bin/metavelvet/meta-velvetg $temp/$reads.velvet -ins_length 380 -min_pair_count 3 -exportFiltered yes -unused_reads yes -alignments yes -min_contig_lgth $contig_length"
done
done
done
done
done
done
done

# DE NOVO ASSEMBLY USING TRINITY
# Trinity is using reads that could not be assembled into contigs with MetaVelvet.
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
for contig_length in 200
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -n 4 -R "span[hosts=1]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.tr "/idi/sabeti-scratch/kandersen/bin/trinity/Trinity.pl --JM 2G --CPU 4 --min_contig_length $contig_length --seqType fa --single $temp/$reads.velvet/UnusedReads.fa --output $temp/$reads.trinity"
done
done
done
done
done
done
done

# CONCATENATE CONTIGS AND REMOVE DUPLICATES & LOW COMPLEXITY
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.rd "cat $temp/$reads.velvet/meta-velvetg.contigs.fa $temp/$reads.trinity/Trinity.fasta | /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -derep 1 -min_gc 10 -max_gc 90 -lc_method dust -lc_threshold 7 -fasta stdin -out_good $url/_meta/$reads.cleaned_contigs -out_bad null"
done
done
done

# PLOT LENGTH OF CONTIGS AND PERFORM VARIOUS SEQUENCE CONTENT ANALYSES - 1
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p1 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -graph_data $temp/$reads.cleaned_contigs.gd -fasta $url/_meta/$reads.cleaned_contigs.fasta -out_good null -out_bad null"
done
done
done

# PLOT LENGTH OF CONTIGS AND PERFORM VARIOUS SEQUENCE CONTENT ANALYSES - 2
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p2 "perl /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-graphs.pl -i $temp/$reads.cleaned_contigs.gd -o $url/_meta/$reads.cleaned_contigs -html_all"
done
done
done

# SPLIT CONTIG FILES
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for num_split_files in 200
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.sc "/idi/sabeti-scratch/kandersen/bin/exonerate/fastasplit -c $num_split_files -f $url/_meta/$reads.cleaned_contigs.fasta -o $temp"
done
done
done
done

# RUN BLASTN ANALYSIS
# Other databases: viral_prokaryot_nt.v4, custom_v1
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt
do
for memory in 10
do
for queue in hour
do
for time in 4:00
do
for evalue in 1e-2
do
for word_size in 16
do
for outformat in 0 # 6
do
i=1
j=1
for a in $temp/$reads.cleaned_contigs.fasta_chunk_*
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bn "blastn -db /idi/sabeti-scratch/kandersen/references/blast/$db -word_size $word_size -evalue $evalue -outfmt $outformat -num_descriptions 50 -num_alignments 50 -query $a -out $temp/$reads.contigs.$db.$((i++)).txt"
done
done
done
done
done
done
done
done
done
done
done

# CONCATENATE BLASTN RESULTS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cb "cat $temp/$reads.contigs.$db.*.txt > $url/_meta/$reads.cleaned_contigs.$db.txt"
done
done
done
done

# CREATE MEGAN FILES
# You will have to be running an xhost in order to do this step: open /Applications/Utilities/X11.app/ xhost + and export DISPLAY=:0. After you have done this, you need to login using ssh -Y $broaduser@cu.broadinstitute.org
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for db in nt
do
for memory in 8
do
for queue in hour
do
for time in 4:00
do
for minSupport in 1
do
for minScore in 0.0
do
for winScore in 0.0
do
for topPercent in 10.0
do
bsub -W 4:00 -R "rusage[mem=8]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.mg /idi/sabeti-scratch/kandersen/bin/megan/MEGAN -E +g -p "/idi/sabeti-scratch/kandersen/bin/megan/Megan.def" -x "load gi2taxfile="/idi/sabeti-scratch/kandersen/bin/megan/class/resources/files/gi_taxid_nucl.bin"; import blastfile="$url/_meta/$reads.cleaned_contigs.$db.txt" fastafile="$url/_meta/$reads.cleaned_contigs.fasta" meganfile="$url/_meta/$reads.cleaned_contigs.$db.rma" useseed=false usekegg=false; recompute minSupport=$minSupport minScore=$minScore topPercent=$topPercent winScore=$winScore minComplexity=0.3 useIdentityFilter=false; quit;"
done
done
done
done
done
done
done
done
done
done

#-------- READ-BASED ANALYSIS --------#
# CONVERT READS TO FASTA AND REMOVE LOW COMPLEXITY
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p1 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -lc_method dust -lc_threshold 7 -out_format 1 -line_width 0 -fastq $url/_reads/$reads.clean.1.fastq -out_good $temp/$reads.cleaned_reads.prinseq.1 -out_bad null"
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p2 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -lc_method dust -lc_threshold 7 -out_format 1 -line_width 0 -fastq $url/_reads/$reads.clean.2.fastq -out_good $temp/$reads.cleaned_reads.prinseq.2 -out_bad null"
done
done
done
done

# FIX MATE-PAIR INFORMATION
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.fm "python /idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastaSeqs.py $temp/$reads.cleaned_reads.prinseq.1.fasta $temp/$reads.cleaned_reads.prinseq.2.fasta $url/_meta/$reads.cleaned_reads"
done
done
done
done

# CONCATENATE READS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -W 4:00 -q hour -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cr "rm $url/_meta/$reads.cleaned_reads.notcommon.fasta && cat $url/_meta/$reads.cleaned_reads.1.fasta $url/_meta/$reads.cleaned_reads.2.fasta > $temp/$reads.clean_concat.fasta"
done
done
done
done

# METAPHLAN ANALYSIS OF BACTERIAL READS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
for analysis_type in rel_ab # reads_map, clade_profiles
do
for evalue in 1e-6
do
for word_size in 16
do
for taxonomic_level in a # k, p, c, o, f, g, s
do
for num_cores in 1
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -n $num_cores -R "span[hosts=1]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.mp "/idi/sabeti-scratch/kandersen/bin/metaphlan/metaphlan.py $temp/$reads.clean_concat.fasta -t $analysis_type --nproc $num_cores --tax_lev $taxonomic_level -o $url/_meta/$reads.cleaned_reads.mphlan.csv --evalue $evalue --word_size $word_size --blastdb /idi/sabeti-scratch/kandersen/bin/metaphlan/blastdb/mpa"
done
done
done
done
done
done
done
done
done
done
done

# SPLIT FILES FOR BLASTN ANALYSIS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.sf "split -a 3 -l 10000 $url/_meta/$reads.cleaned_reads.1.fasta $temp/$reads.cleaned_reads.1.split."
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.sf "split -a 3 -l 10000 $url/_meta/$reads.cleaned_reads.2.fasta $temp/$reads.cleaned_reads.2.split."
done
done
done

# RUN BLASTN ANALYSIS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt # custom_v1, viral_prokaryot_nt.v4
do
for evalue in 1e-2
do
for word_size in 16
do
for outformat in 0 # 6
do
for memory in 10
do
i=1
j=1
for a in $temp/$reads.cleaned_reads.1.split.*
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bn "blastn -db /idi/sabeti-scratch/kandersen/references/blast/$db -word_size $word_size -evalue $evalue -outfmt $outformat -num_descriptions 50 -num_alignments 50 -query $a -out $temp/$reads.1.$db.$((i++)).txt"
done
done
done
done
done
done
done
done
done
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt # custom_v1, viral_prokaryot_nt.v4
do
for evalue in 1e-2
do
for word_size in 16
do
for outformat in 0 # 6
do
for memory in 10
do
i=1
j=1
for b in $temp/$reads.cleaned_reads.2.split.*
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bn "blastn -db /idi/sabeti-scratch/kandersen/references/blast/$db -word_size $word_size -evalue $evalue -outfmt $outformat -num_descriptions 50 -num_alignments 50 -query $b -out $temp/$reads.2.$db.$((i++)).txt"
done
done
done
done
done
done
done
done
done

# CONCATENATE BLASTN RESULTS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt # custom_v1, viral_prokaryot_nt.v4
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cb "cat $temp/$reads.?.$db.*.txt > $url/_meta/$reads.cleaned_reads.$db.txt"
done
done
done
done

# CREATE MEGAN FILES
# You will have to be running an xhost in order to do this step: open /Applications/Utilities/X11.app/ xhost + and export DISPLAY=:0. After you have done this, you need to login using ssh -Y $broaduser@cu.broadinstitute.org
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt # custom_v1, viral_prokaryot_nt.v4
do
for memory in 8
do
for queue in hour
do
for time in 4:00
do
for minSupport in 1
do
for minScore in 0.0
do
for winScore in 0.0
do
for topPercent in 10.0
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.mg /idi/sabeti-scratch/kandersen/bin/megan/MEGAN -E +g -p "/idi/sabeti-scratch/kandersen/bin/megan/Megan.def" -x "load gi2taxfile="/idi/sabeti-scratch/kandersen/bin/megan/class/resources/files/gi_taxid_nucl.bin"; import blastfile="$url/_meta/$reads.cleaned_reads.$db.txt" fastafile='$url/_meta/$reads.cleaned_reads.1.fasta', '$url/_meta/$reads.cleaned_reads.2.fasta' meganfile="$url/_meta/$reads.cleaned_reads.$db.rma" useseed=false usekegg=false paired=true suffix1='/1' suffix2='/2'; recompute minSupport=$minSupport minScore=$minScore topPercent=$topPercent winScore=$winScore minComplexity=0.3 useIdentityFilter=false; quit;"
done
done
done
done
done
done
done
done
done
done
done

# EXTRACT READS WITH NO BLAST HITS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for db in nt # custom_v1, viral_prokaryot_nt.v4
do
for script in noBlastHits_v2.py # use noBlastHits.py if you used tabular (6) format in blast search
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.nb "/idi/sabeti-scratch/kandersen/bin/scripts/$script -b $url/_meta/$reads.cleaned_reads.$db.txt -r $url/_meta/$reads.cleaned_reads.1.fasta -m nohit > $temp/$reads.nohits.1.fasta"
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.nb "/idi/sabeti-scratch/kandersen/bin/scripts/$script -b $url/_meta/$reads.cleaned_reads.$db.txt -r $url/_meta/$reads.cleaned_reads.2.fasta -m nohit > $temp/$reads.nohits.2.fasta"
done
done
done
done
done

# FIX MATE-PAIR INFORMATION
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.fm "python /idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastaSeqs.py $temp/$reads.nohits.1.fasta $temp/$reads.nohits.2.fasta $url/_meta/$reads.nohits_reads"
done
done
done
done

#-------- ANALYZE NO-HIT READS --------#
# CONCATENATE NO-HIT READS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
bsub -W 4:00 -q hour -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cr "rm $url/_meta/$reads.nohits_reads.notcommon.fasta && cat $url/_meta/$reads.nohits_reads.?.fasta > $temp/$reads.nohits.cat.fasta"
done
done
done
done

# PCA AND OTHER ANALYSES
input: $url/_meta/$reads.nohits_reads.1.fasta / $url/_meta/$reads.nohits_reads.2.fasta

# SPLIT READS FOR BLASTX ANALYSIS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.sf "split -a 3 -l 1000 $temp/$reads.nohits.cat.fasta $temp/$reads.9a.split.fasta."
done
done
done

# RUN BLASTX ON NOVEL READS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for queue in hour
do
for time in 4:00
do
for evalue in 1e2
do
for db in nr
do
i=1
j=1
for a in $temp/$reads.9a.split.fasta.*
do
bsub -R "rusage[mem=8]" -q $queue -W $time -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bx "blastx -db /idi/sabeti-scratch/kandersen/references/blast/$db -evalue $evalue -outfmt 0 -num_descriptions 50 -num_alignments 50 -query $a -out $temp/$reads.15a.blastx.$((i++)).txt"
done
done
done
done
done
done
done
done

# CONCATENATE BLASTX RESULTS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cx "cat $temp/$reads.15a.blastx.*.txt > $url/_meta/$reads.nohits_reads.nr.txt"
done
done
done

# CREATE MEGAN FILES
# You will have to be running an xhost in order to do this step: open /Applications/Utilities/X11.app/ xhost + and export DISPLAY=:0. After you have done this, you need to login using ssh -Y $broaduser@cu.broadinstitute.org
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 8
do
for queue in hour
do
for time in 4:00
do
for minSupport in 1
do
for minScore in 0.0
do
for winScore in 0.0
do
for topPercent in 10.0
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.mg /idi/sabeti-scratch/kandersen/bin/megan/MEGAN -E +g -p "/idi/sabeti-scratch/kandersen/bin/megan/Megan.def" -x "load gi2taxfile="/idi/sabeti-scratch/kandersen/bin/megan/class/resources/files/gi_taxid_nucl.bin"; import blastfile="$url/_meta/$reads.nohits_reads.nr.txt" fastafile="$temp/$reads.nohits.cat.fasta" meganfile="$url/_meta/$reads.nohits_reads.nr.rma" useseed=true usekegg=true; recompute minSupport=$minSupport minScore=$minScore topPercent=$topPercent winScore=$winScore minComplexity=0.3 useIdentityFilter=false; quit;"
done
done
done
done
done
done
done
done
done
done

#-------- CONTIG ANALYSIS OF NO-HIT READS --------#
# DE NOVO ASSEMBLY USING TRINITY
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
for contig_length in 200
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -n 4 -R "span[hosts=1]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.dn "/idi/sabeti-scratch/kandersen/bin/trinity/Trinity.pl --JM 2G --CPU 4 --min_contig_length $contig_length --seqType fa --left $url/_meta/$reads.nohits_reads.1.fasta --right $url/_meta/$reads.nohits_reads.2.fasta --output $temp/$reads.trinity_novel"
done
done
done
done
done
done
done

# RENAME AND COPY CONTIGS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cc "mv $temp/$reads.trinity_novel/Trinity.fasta $url/_meta/$reads.nohits_contigs.fasta"
done
done
done

#-------- RUN ANALYSIS ON NOVEL CONTIGS - ONLY DO THIS IF YOU FOUND CONTIGS --------#
# PLOT LENGTH OF CONTIGS AND PERFORM VARIOUS SEQUENCE CONTENT ANALYSES - 1
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p1 "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -graph_data $temp/$reads.cleaned_reads.nohits.contigs.gd -fasta $url/_meta/$reads.nohits_contigs.fasta -out_good null -out_bad null"
done
done
done

# SPLIT FASTA FILE FOR BLASTX ANALYSIS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for num_split in 50
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.sx "/idi/sabeti-scratch/kandersen/bin/exonerate/fastasplit -c $num_split -f $url/_meta/$reads.nohits_contigs.fasta -o $temp/"
done
done
done
done

# PLOT LENGTH OF CONTIGS AND PERFORM VARIOUS SEQUENCE CONTENT ANALYSES - 2
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.p2 "perl /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-graphs.pl -i $temp/$reads.cleaned_reads.nohits.contigs.gd -o $url/_meta/$reads.nohits_contigs -html_all"
done
done
done

# RUN BLASTX ON NOVEL CONTIGS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
for queue in hour
do
for time in 4:00
do
for evalue in 1e2
do
for db in nr
do
i=1
j=1
for a in $temp/$reads.nohits_contigs.fasta_chunk_*
do
bsub -R "rusage[mem=8]" -q $queue -W $time -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.$((j++)).bx "blastx -db /idi/sabeti-scratch/kandersen/references/blast/$db -evalue $evalue -outfmt 0 -num_descriptions 50 -num_alignments 50 -query $a -out $temp/$reads.15.blastx.$((i++)).txt"
done
done
done
done
done
done
done
done

# CONCATENATE BLASTX RESULTS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cx "cat $temp/$reads.15.blastx.*.txt > $url/_meta/$reads.nohits_contigs.nr.txt"
done
done
done

# CREATE MEGAN FILES
# You will have to be running an xhost in order to do this step: open /Applications/Utilities/X11.app/ xhost + and export DISPLAY=:0. After you have done this, you need to login using ssh -Y $broaduser@cu.broadinstitute.org
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for memory in 8
do
for queue in hour
do
for time in 4:00
do
for minSupport in 1
do
for minScore in 0.0
do
for winScore in 0.0
do
for topPercent in 10.0
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_meta -G sabeti_labfolk -J $reads.mg /idi/sabeti-scratch/kandersen/bin/megan/MEGAN -E +g -p "/idi/sabeti-scratch/kandersen/bin/megan/Megan.def" -x "load gi2taxfile="/idi/sabeti-scratch/kandersen/bin/megan/class/resources/files/gi_taxid_nucl.bin"; import blastfile="$url/_meta/$reads.nohits_contigs.nr.txt" fastafile="$url/_meta/$reads.nohits_contigs.fasta" meganfile="$url/_meta/$reads.nohits_contigs.nr.rma" useseed=true usekegg=true; recompute minSupport=$minSupport minScore=$minScore topPercent=$topPercent winScore=$winScore minComplexity=0.3 useIdentityFilter=false; quit;"
done
done
done
done
done
done
done
done
done

#-------- FINISH, CLEAN AND MOVE RELEVANT FILES --------#
# CALCULATE NUMBER OF READS IN EACH STEP
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/$reads.log.counts.txt -P sabeti_meta -G sabeti_labfolk -J $reads.ca "echo 'Order is: Start:Quality:BT_transcript:BT_hg19:BT_Contaminants:Stampy:Duplicates:Complexity:NoBlastHits:Contigs:NoContigs:NoHitContigs' && cat $url/_reads/$reads.reads?.fastq | wc -l && cat $temp/$reads.trimmed.?.fastq | wc -l && cat $temp/$reads.bmtagger.mrna.?.fastq | wc -l && cat $temp/$reads.bmtagger.hg19.?.fastq | wc -l && cat $temp/$reads.bmtagger.contaminants.?.fastq | wc -l && cat $url/_reads/$reads.clean.?.fastq | wc -l && cat $url/_meta/$reads.cleaned_reads.?.fasta | wc -l && cat $url/_meta/$reads.nohits_reads.?.fasta | wc -l && cat $url/_meta/$reads.cleaned_contigs.fasta | egrep -e '>' | wc -l && cat $temp/$reads.velvet/UnusedReads.fa | egrep -e '>' | wc -l && cat $url/_meta/$reads.nohits_contigs.fasta | egrep -e '>' | wc -l"
done
done
done

# COMPRESS CLEANED READS & CONTIGS
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
bsub -o $url/_logs/0cleanup.txt -P sabeti_meta -G sabeti_labfolk -J $reads.cs "gzip -9 $url/_meta/$reads.*.fasta"
done
done

# CLEANUP
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for temp in /broad/hptmp/stremlau
do
bsub -o $url/_logs/0cleanup.txt -P sabeti_meta -G sabeti_labfolk -J $reads.rm "rm -rf $temp/$reads* && rm -rf $url/_reads/$reads* && rm -rf $url/_meta/$reads.*.txt && rm -rf $url/_reports/$reads.*_fastqc"
done
done
done

# MOVE RELEVANT FILES
for reads in 49CMiSeq
do
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for target in /idi/sabeti-data/kandersen/metagenomic_storage
do
bsub -W 4:00 -o $url/_logs/0cleanup.txt -P sabeti_align -G sabeti_labfolk -J $reads.mv "mv $url/_meta/$reads.*.fasta.gz $target/cleaned_fasta/ && mv $url/_meta/$reads.*.rma $target/megan/ && mv $url/_reports/$reads.xlsx $target/reports_overall/ && mv $url/_meta/$reads.*.html $target/reports_prinseq/ && mv $url/_reports/$reads.*.zip $target/reports_quality/"
done
done
done