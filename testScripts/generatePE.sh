#!/bin/bash
# This script converts a series of mRNA sequencing data file in FASTQ format
# to a table of UMI read counts of human genes in multiple sample conditions.

# 1 Parameters

# 1.1 Global

TOP_DIR=$1

# 1.2 Dataset
SERIES="20150409"
SAMPLE_ID="RNAseq_${SERIES}"
LANES=6
DATA_DIR=$TOP_DIR/data/LINCS
SEQ_DIR="${DATA_DIR}/peSeqs"
ALIGN_DIR="${DATA_DIR}/Aligns"
COUNT_DIR="${DATA_DIR}/Counts_sam"
COUNT_DIR_SAF="${DATA_DIR}/Counts_saf"

mkdir -p $COUNT_DIR
mkdir -p $COUNT_DIR_SAF

UMITOOLS_DIR="${TOP_DIR}"
UMITOOLS_OLD_DIR="/home/lhhung/LINCS_RNAseq_cpp"
REF_DIR="$DATA_DIR/References/Broad_UMI"
SPECIES_DIR="${REF_DIR}/Human_RefSeq"
REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa"
BARCODE_FILE="${REF_DIR}/barcodes_trugrade_96_set4.dat"
FASTQ_LIST="${DATA_DIR}/fastqList"
filterCmd="${UMITOOLS_DIR}/source/w96/umimerge_filter -B 0 -s $SYM2REF_FILE -e $ERCC_SEQ_FILE -S -U 0"
# 1.4 Program
PROG_DIR="$DATA_DIR/Programs/Broad-DGE"
BWA_ALN_SEED_LENGTH=24
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=20
THREAD_NUMBER=8

# 2 Computation

# 2.1 Alignment
# Align sequence fragments to reference genome library.

 SEQ_FILES="${SEQ_DIR}/SRR6799772_1.fastq ${SEQ_DIR}/SRR6799772_2.fastq"

 #use tight checking no mismatch no ambiguities to match original - default is the looser setting of mismatch =1 and missing N=1 

echo "$UMITOOLS_DIR/source/w96/split_pe -v -s 250 -o $ALIGN_DIR -t $THREAD_NUMBER  $SEQ_FILES"
$UMITOOLS_DIR/source/w96/split_pe -v -s 25000 -o $ALIGN_DIR -t $THREAD_NUMBER  $SEQ_FILES
echo "${UMITOOLS_DIR}/testScripts/multibwaPE.sh $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER $filterCmd"
$UMITOOLS_DIR/testScripts/multibwaPE.sh $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER "\"${filterCmd}\""
exit
#echo "${UMITOOLS_DIR}/testScripts/multibwaNoMulti.sh $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER"
#$UMITOOLS_DIR/testScripts/multibwaNoMulti.sh $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER 


echo "$UMITOOLS_DIR/source/w96/umimerge_parallel -i $SAMPLE_ID -b $FASTQ_LIST -e $ERCC_SEQ_FILE -s $SYM2REF_FILE -a $ALIGN_DIR -o $COUNT_DIR -t $THREAD_NUMBER -F 6 -S -U 0" 
s$UMITOOLS_DIR/source/w96/umimerge_parallel -i $SAMPLE_ID -b $FASTQ_LIST -e $ERCC_SEQ_FILE -s $SYM2REF_FILE -a $ALIGN_DIR -o $COUNT_DIR -t $THREAD_NUMBER -F 6 -S -U 0 

echo "$UMITOOLS_DIR/source/w96/umimerge_parallel -i $SAMPLE_ID -b $FASTQ_LIST -e $ERCC_SEQ_FILE -s $SYM2REF_FILE -a $ALIGN_DIR -o $COUNT_DIR_SAF -t $THREAD_NUMBER -f -F 6 -S -U 0" 
$UMITOOLS_DIR/source/w96/umimerge_parallel -i $SAMPLE_ID -b $FASTQ_LIST -e $ERCC_SEQ_FILE -s $SYM2REF_FILE -a $ALIGN_DIR -o $COUNT_DIR_SAF -t $THREAD_NUMBER -f -F 6 -S -U 0 
