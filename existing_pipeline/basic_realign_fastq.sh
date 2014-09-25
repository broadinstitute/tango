1) download reference sequence
2) put into Geneious and then convert to fasta
3) put it into reference folder and run following script

for species in 
do
for reference in 
do
bsub -o /idi/sabeti-scratch/kandersen/references/matt/rhabdo/index.log.bsub.txt -P sabeti_kga -J $reference.index "/idi/sabeti-scratch/kandersen/bin/novocraft/novoindex /idi/sabeti-scratch/kandersen/references/$species/$reference.nix /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta && bwa index -p /idi/sabeti-scratch/kandersen/references/$species/$reference /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar REFERENCE=/idi/sabeti-scratch/kandersen/references/$species/$reference.fasta O=/idi/sabeti-scratch/kandersen/references/$species/$reference.dict"
done
done


#-------- GENOME RE-ALIGNMENT @ BROAD --------#
# NO TRAILING / IN URL
use Perl-5.10
use Samtools
use Python-2.7
use Java-1.6

# RE-ALIGN THE READS
# coding: ILMFQ (1.5); STDFQ (1.8); FA (fasta). Date: yymmdd, e.g. 120124. SeqCenter: Harvard/Broad.
# For less strict alignment use -l 25 -g 40 -x 20 or completely leave out
for reads in
do
for url in
do
for species in matt/
do
for reference in 
do
for coding in STDFQ
do
for date in
do
for library in BroadPE
do
for seq_technology in Illumina
do
for platform in HiSeq2k
do
for SeqCenter in Broad
do
for suffix in fastq
do
for memory in 2
do
for queue in hour
do
for time in 4:00
do
bsub -W $time -q $queue -R "rusage[mem=$memory]" -n 4 -R "span[hosts=1]" -o $url/_logs/$reads.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J $reads.al "/idi/sabeti-scratch/kandersen/bin/novocraft/novoalign -k -c 2 -f $url/_reads/$reads.reads1.$suffix $url/_reads/$reads.reads2.$suffix -r Random -t 500 -F $coding -d /idi/sabeti-scratch/kandersen/references/$species/$reference.nix -o SAM $'@RG\tID:$date.$reads\tSM:$reads\tPL:$seq_technology\tPU:$platform\tLB:$library\tCN:$SeqCenter' 2> $url/_logs/$reads.log.novoalign.txt | samtools view -b -S -q 1 -u - | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$url/_bams/$reads.mapped.bam CREATE_INDEX=true && java -Xmx2g -jar /seq/software/picard/current/bin/MarkDuplicates.jar I=$url/_bams/$reads.mapped.bam O=$url/_bams/$reads.mappedNoDub.bam METRICS_FILE=$url/_logs/$reads.log.markdups.txt CREATE_INDEX=true REMOVE_DUPLICATES=true"
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
done
done
done

# CALCULATE COVERAGE, CONSENSUS AND SNPS
for reads in
do
for url in
do
for species in lassa
do
for reference in $reads
do
for suffix in mapped.bam
do
bsub -W 4:00 -o $url/_logs/$reads.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J $reads.st "/idi/sabeti-scratch/kandersen/bin/bedtools/genomeCoverageBed -d -ibam $url/_bams/$reads.$suffix -g /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta > $url/_pileup/$reads.coverage.txt && samtools mpileup -Q 15 -uB -q 1 -d 10000 -f /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta $url/_bams/$reads.$suffix | bcftools view -Acg - | /idi/sabeti-scratch/kandersen/bin/scripts/vcfutils.pl vcf2fq > $url/_pileup/$reads.cns.fastq && samtools mpileup -Q 15 -B -q 1 -d 10000 -f /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta $url/_bams/$reads.$suffix | java -jar /idi/sabeti-scratch/kandersen/bin/varscan/varscan.jar pileup2snp --min-reads2 4 > $url/_pileup/$reads.snps.txt && bamtools stats -insert -in $url/_bams/$reads.$suffix > $url/_logs/$reads.log.bamstats_realigned.txt"
done
done
done
done
done