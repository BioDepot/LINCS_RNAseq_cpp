#lhhung 0319
#!/bin/bash
# ARGS are $TOP_DIR,$REF_DIR,$SPECIES_DIR,$ALIGN_DIR,$BWA_ALN_SEED_LENGTH,$BWA_SAM_MAX_ALIGNS_FOR_XA_TAG,FILTERBIN,FILTERFLAGS,nThreads
TOP_DIR=$1
REF_DIR=$2
SPECIES_DIR=$3
ALIGN_DIR=$4
BWA_ALN_SEED_LENGTH=$5
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=$6
nThreads=$7
FILTERBIN=${8:1:-1}
suffix="saf"
 if [ -z "$FILTERBIN" ]; then
  FILTERBIN="grep -v '^\@'"
  suffix="sam"
 fi

REF_SEQ_FILE=$SPECIES_DIR/refMrna_ERCC_polyAstrip.hg19.fa


lockDir=/tmp/locks.$$
mkdir -p $lockDir

runJob(){
	#pid=$( sh -c 'echo $PPID' )
	lasti=$((${#files[@]} - 1))
	#lasti=0
 for i in $(seq 0 ${lasti}); do
  if (mkdir $lockDir/lock$i 2> /dev/null ); then
   f1=${files[$i]}
   f1begin=${f1%1_*}
   f1end=${f1##*1_}
   f2="$f1begin"2_"$f1end"
   echo thread $1 working on $f1 $f2
   echo "(bwa aln -l $BWA_ALN_SEED_LENGTH -t 1 $REF_SEQ_FILE $f1 $f2 2>/dev/null | bwa samse -n $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $REF_SEQ_FILE - $file | $FILTERBIN > $file.$suffix )" 		  
   (bwa aln -l $BWA_ALN_SEED_LENGTH -t 1 $REF_SEQ_FILE $f1 $f2 2>/dev/null | bwa samse -n $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $REF_SEQ_FILE - $file | $FILTERBIN > $f1.$suffix ) 
   (bwa aln -l $BWA_ALN_SEED_LENGTH -t 1 $REF_SEQ_FILE $f1 $f2 2>/dev/null | bwa samse -n $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $REF_SEQ_FILE - $file | grep -v '^\@' > $f1.sam)
   #rm $file
 fi
done
exit
}

files=( $(find $ALIGN_DIR -mindepth 1 -maxdepth 2 -name "*1_0*.fq" -type f ) )

for i in $(seq 2 $nThreads); do
	runJob $i &
done

runJob 1 &
wait
rm -rf $lockDir
