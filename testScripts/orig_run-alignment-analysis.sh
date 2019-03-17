#!/bin/bash
# This script converts a series of mRNA sequencing data file in FASTQ format
# to a table of UMI read counts of human genes in multiple sample conditions.

#modified to use cpp program Hong Hung 2016

# 1.1 Global
TOP_DIR="/storage/detoxs"

# 1.2 Dataset
SERIES="20150409"
SAMPLE_ID="RNAseq_${SERIES}"
LANES=6
DATA_DIR=$TOP_DIR/UMITest/
SEQ_DIR="${DATA_DIR}/Seqs"
ALIGN_DIR="${DATA_DIR}/Aligns"
COUNT_DIR="${DATA_DIR}/Counts"
UMITOOLS_DIR="$TOP_DIR/LINCS_RNAseq_cpp/" #where the cpp binaries are
# 1.3 Reference
REF_DIR="$TOP_DIR/Broad_UMI"
SPECIES_DIR="${REF_DIR}/Human_RefSeq"
REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa"
BARCODE_FILE="${REF_DIR}/barcodes_trugrade_96_set4.dat"

# 1.4 Program
PROG_DIR="$TOP_DIR/Programs/Broad-DGE"
BWA_ALN_SEED_LENGTH=24
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=20
THREAD_NUMBER=6

# 2 Computation

# 2.1 Alignment
# Align sequence fragments to reference genome library.
let "IDX = 1"
while [ "$IDX" -le "${LANES}" ]; do
	SUBSAMPLE_ID="Lane$IDX"
	SEQ_FILE_R1="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq"
	SEQ_FILE_R2="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq"
	SEQ_FILES="${SEQ_FILES} ${SEQ_FILE_R1} ${SEQ_FILE_R2}"
	let "IDX = $IDX + 1"
done	

#echo "${UMITOOLS_DIR}/umisplit_orig -v -l 16 -t $THREAD_NUMBER -o $ALIGN_DIR $SEQ_FILES"
#$UMITOOLS_DIR/umisplit_orig -v -l 16 -t $THREAD_NUMBER -o $ALIGN_DIR $SEQ_FILES

#echo "${UMITOOLS_DIR}/scripts/multibwa_orig.pl $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER"
#$UMITOOLS_DIR/scripts/multibwa_orig.pl $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER
echo "$UMITOOLS_DIR/umimerge $SAMPLE_ID $SYM2REF_FILE $ERCC_SEQ_FILE $BARCODE_FILE $ALIGN_DIR $COUNT_DIR FALSE"
$UMITOOLS_DIR/umimerge $SAMPLE_ID $SYM2REF_FILE $ERCC_SEQ_FILE $BARCODE_FILE $ALIGN_DIR $COUNT_DIR FALSE

exit 0
