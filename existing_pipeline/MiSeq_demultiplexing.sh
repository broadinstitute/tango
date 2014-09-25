#-------- DEMULTIPLEXING OF MISEQ AND HISEQ DATA @ BROAD --------#
# Picard directory e.g. /seq/picard/D0N2CACXX/C1-210_2012-02-29_2012-04-12/2
# Bustard directory e.g. /seq/solexaproc/SL-HAW/analyzed/120229_SL-HAW_0168_AFCD0N2CACXX
# Bustard directory can be found by looking in the following file: /seq/picard/<flowcell>/<analysis>/info/logs/ID_xxx.json - look for 'RunFolder'
# Barcode for ExtractIlluminaBarcodes e.g. /seq/solexaproc/SL-HAW/analyzed/120229_SL-HAW_0168_AFCD0N2CACXX/Data/Intensities/BaseCalls/barcodeData.2
# Barcode for IlluminaBasecallsToSam e.g. /seq/picard/D0N2CACXX/C1-210_2012-02-29_2012-04-12/2/library_params.txt
# ExtractIlluminaBarcodes barcode example: http://cl.ly/3g1M2D0C2b2C3C0M0h2C
# IlluminaBasecallsToSam barcode example: http://cl.ly/1N2w2n2U3Z431D382Q3j

# MAKE REQUIRED SUB-DIRECTORIES - DON'T DO THIS IF YOU ALREADY CREATED THESE WITH ANOTHER PIPELINE
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
bsub -o ~/log.txt -P sabeti_align -G sabeti_labfolk -J Directories "mkdir $url/_logs $url/_temp $url/_reads $url/_bams $url/_reports $url/_pileup $url/_meta"
done

# CREATE LINKS TO RAW DATA
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for lane in 1
do
for date in 130110
do
for flow_cell in 000000000-A2KP4
do
for temp in /broad/hptmp/stremlau
do
for bustard_dir in /idi/sabeti-data/stremlau/130110_49CMiSeq/Raw/130107_M00266_0125_000000000-A2KP4
do
bsub -q hour -o $url/_logs/barcodes.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J bc.sl "mkdir $temp/$date.$flow_cell && mkdir $temp/$date.$flow_cell/Data && mkdir $temp/$date.$flow_cell/Data/Intensities && mkdir $temp/$date.$flow_cell/Data/Intensities/BaseCalls && ln -s $bustard_dir/Data/Intensities/L00* $temp/$date.$flow_cell/Data/Intensities && ln -s $bustard_dir/Data/Intensities/BaseCalls/L00* $temp/$date.$flow_cell/Data/Intensities/BaseCalls && ln -s $bustard_dir/Data/Intensities/BaseCalls/barcodeData* $temp/$date.$flow_cell/Data/Intensities/BaseCalls"
done
done
done
done
done
done

# EXTRACT ILLUMINA BARCODES
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq 
do
for lane in 1
do
for date in 130110
do
for flow_cell in 000000000-A2KP4
do
for temp in /broad/hptmp/stremlau
do
for read_structure in 301T8B151T # 151T8B151T 151T6B151T
do
for max_mismatches in 0
do
for min_quality_score in 25
do
for bar_codes in $temp/$date.$flow_cell/Data/Intensities/BaseCalls/barcodeData.$lane
do
for memory in 4
do
bsub -n 4,8 -R "span[hosts=1] rusage[mem=$memory]" -q hour -W 4:00 -o $url/_logs/barcodes.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J bc.eb "java -Xmx2g -jar /seq/software/picard/current/bin/ExtractIlluminaBarcodes.jar BASECALLS_DIR=$temp/$date.$flow_cell/Data/Intensities/BaseCalls/ LANE=$lane READ_STRUCTURE=$read_structure BARCODE_FILE=$bar_codes METRICS_FILE=$url/_logs/barcode.metrics.txt MAX_MISMATCHES=$max_mismatches MINIMUM_BASE_QUALITY=$min_quality_score NUM_PROCESSORS=0"
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

# CONVERT RAW FILES TO BAM
for url in /idi/sabeti-data/stremlau/130110_49CMiSeq
do
for lane in 1
do
for date in 130110
do
for flow_cell in 000000000-A2KP4
do
for temp in /broad/hptmp/stremlau
do
for read_structure in 101T8B101T # 151T8B151T 151T6B151T
do
for library_parameters in $url/_logs/library_params.txt
do
bsub -q hour -W 4:00 -n 4,8 -R "span[hosts=1] rusage[mem=16]" -o $url/_logs/barcodes.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J bc.cs "java -Xmx12g -jar /seq/software/picard/current/bin/IlluminaBasecallsToSam.jar BASECALLS_DIR=$temp/$date.$flow_cell/Data/Intensities/BaseCalls/ LANE=$lane READ_STRUCTURE=$read_structure LIBRARY_PARAMS=$library_parameters SEQUENCING_CENTER=Broad RUN_BARCODE=$flow_cell NUM_PROCESSORS=0 ADAPTERS_TO_CHECK=PAIRED_END MAX_READS_IN_RAM_PER_TILE=100000 MAX_RECORDS_IN_RAM=300000 FORCE_GC=false"
done
done
done
done
done
done
done

# CLEANUP
for date in 130110
do
for flow_cell in 000000000-A2KP4
do
for temp in /broad/hptmp/stremlau
do
bsub -q hour -W 4:00 -o $url/_logs/barcodes.log.bsub.txt -P sabeti_align -G sabeti_labfolk -J bc.cl "rm -rf $temp/$date.$flow_cell"
done
done
done