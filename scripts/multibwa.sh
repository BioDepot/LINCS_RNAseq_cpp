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
FILTERBIN="grep -v '^\@'"
if [ -z "$8" ]; then
 FILTERBIN="grep -v '^\@'"
else
 eval FILTERBIN = $8
fi
REF_SEQ_FILE=$SPECIES_DIR/refMrna_ERCC_polyAstrip.hg19.fa

lockDir=/tmp/locks.$$
mkdir -p $lockDir

runJob(){
	#pid=$( sh -c 'echo $PPID' )
	lasti=$((${#dirs[@]} - 1))
 for i in $(seq 0 ${lasti}); do
  if (mkdir $lockDir/lock$i 2> /dev/null ); then
   dir=${dirs[$i]}
		 echo thread $1 working on $dir
		 #echo "cat $dir/*.fq >  $dir/all.fastq"
		 cat $dir/*.fq >  $dir/all.fastq
		 (bwa aln -l $BWA_ALN_SEED_LENGTH -t 1 $REF_SEQ_FILE $dir/all.fastq 2>/dev/null | bwa samse -n $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $REF_SEQ_FILE - $dir/all.fastq | $FILTERBIN > $dir/all.sam ) & > /dev/null 
		 rm $dir/all.fastq
		fi
	done
	exit
}
dirs=( $(find $ALIGN_DIR -mindepth 1 -maxdepth 1 -type d))

for i in $(seq 2 $nThreads); do
	  runJob $i &
done
runJob 1 &
wait
rm -rf $lockDir
