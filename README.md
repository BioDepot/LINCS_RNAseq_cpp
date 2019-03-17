
# Cpp based tools for barcoded RNAseq

This a set of cpp binaries and scripts for (Unique Molecular Identifier) UMI based RNAseq. The pre-print is available [here.](https://www.biorxiv.org/content/10.1101/345819v2)  These are decendants of python scripts describe  [here.](https://www.biorxiv.org/content/early/2014/03/05/003236)  (https://www.biorxiv.org/content/early/2014/03/05/003236) used to split fastq files before alignment and merge files after alignment by bwa. All the binaries are multithreaded unless otherwise stated.

## To compile and run
### Clone repo
	git clone https://github.com
### Compile executable
	cd LINCS_RNAseq_cpp/source
	make clean; make all96; make all384
This will result in executable binaries created in the LINCS_RNAseq_cpp/source/w96 (for 96 plate wells) and LINCS_RNAseq_cpp/source/w384 directories (for 384 plate wells).
### Run example
	cd LINCS_RNAseq_cpp
	./runExample.sh
	
This script download all the necessary support files (reference indices and barcode files) and a small dataset. The transcript counts will be in the LINCS_RNAseq_cpp/data/LINCS/Counts directory

For all binaries the -h flag gives documentation about the available flags and an example of how to use the binary: 

##To run using Docker

As an alternative to compiling and running the executables, we also provide a docker container with the executables installed.

### Building Docker containers locally

As an alternative to pulling the Docker container from our Biodepot repository, the user can also build the containers locally using the provided Docker files. 

	cd  LINCS_RNAseq_cpp
	scripts/build.sh

If the container is not built, it will be automatically downloaded when the demo script is run.

### Running the example using Docker

	cd LINCS_RNAseq_cpp
	./runDockerExample.sh

## How to adapt the software for user data

### Using local executables 
The executables themselves can be used in any way that the user wishes. See the documentation for the individual binaries. Binaries for 96 plate wells are created in the w96 directory, whereas the w384 directory contains executables when using 384 wells.
On a Fedora distro, on the example data and running from the github directory - the command for umisplit on a 96 well plate would be

	./source/w96/umisplit -v -l 16 -m 0 -N 0 -o ./data/LINCS/Aligns  -b ./data/LINCS/References/Broad_UMI/barcodes_trugrade_96_set4.dat -t 4


### Using the Docker container
The Docker container can be used if the user is not in a Linux environment and does not wish to compile the the cpp code.

Docker must be installed and running. Then the Docker container is invoked before the command. The local directory needs to be mapped to an internal docker directory and an environment variable must be passed indicating the number of wells.  
On a Fedora distro - the command for umisplit on a 96 well plate would be:

	sudo service docker start
	docker run --rm -v ${PWD}/data/LINCS:/data -e NWELLS=96 umisplit -v -l 16 -m 0 -N 0 -o /data/Aligns  -b /data/References/Broad_UMI/barcodes_trugrade_96_set4.dat -t 4

### Customizing the shell scripts as a guide for user data 

Customizing the parameters in theshell scripts provided for local execution and execution with Docker is probably the easiest method to apply the software to the user data

For local execution the script to call is fast_run_alignment_analysis.sh in the scripts directory of the GitHub repo

To use, the command to type (from the gitHub directory is)

	scripts/fast_run_alignment_analysis.sh ${PWD}
	
To customize the script - open up the script in an editor and modify the parameters


	SERIES="20150409"  <-- used to construct the filenames in the script
	SAMPLE_ID="RNAseq_${SERIES}" <-- used to construct the filenames
	LANES=6 <-- used in loop to construct sample filenames
	DATA_DIR=$TOP_DIR/data/LINCS <-- where the data resides relative (the TOP_DIR variable is where the script is being invoked from)
	SEQ_DIR="${DATA_DIR}/Seqs" <-- where the sequences reside
	ALIGN_DIR="${DATA_DIR}/Aligns" <-- where the alignment is being output
	COUNT_DIR="${DATA_DIR}/Counts" <-- where the Counts will be generated (the Counts/<filename>.unq.refseq.umi.dat file contains the UMI filtered transcript counts) 
	UMITOOLS_DIR="${TOP_DIR}"
	REF_DIR="$DATA_DIR/References/Broad_UMI" <-- the GitHub repo directory - used to calculate where the executables reside
	SPECIES_DIR="${REF_DIR}/Human_RefSeq" <-- This is where the reference indices, translation tables for transcript and spike-ins are stored
	REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa" <-- reference sequence (data is actually unused by bwa but necessary to calculate the file paths for the index files)
	SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat" <-- translation file for transcript names to gene names
	ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa" <-- spike in file 
	BARCODE_FILE="${REF_DIR}/barcodes_trugrade_96_set4.dat" <-- barcodes file
	
In addition to changing these parameters. The following loop can be customized to construct the user fastq names

	let "IDX = 1"	
	SEQ_FILES="";
	while [ "$IDX" -le "${LANES}" ]; do
		SUBSAMPLE_ID="Lane$IDX"
		SEQ_FILE_R1="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq.gz"
		SEQ_FILE_R2="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq.gz"
		SEQ_FILES="${SEQ_FILES} ${SEQ_FILE_R1} ${SEQ_FILE_R2}"
		let "IDX = $IDX + 1"
	done	

Finally the executable commands in the script assume 96 well plates. The commands w96 needs to be changed to w384 if 384 well plates are being used i.e.

	$UMITOOLS_DIR/source/w96/umisplit -v -l 16 -m 0 -N 0 -f -o $ALIGN_DIR -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES 
	
would be changed to

	$UMITOOLS_DIR/source/w384/umisplit -v -l 16 -m 0 -N 0 -f -o $ALIGN_DIR -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES 

The Docker script is invoked by

	scripts/docker_fast_run_alignment_analysis.sh ${PWD}

It is structured in the same way. However, all paths are relative to the internal container file paths. 

The shell scripts are meant as a guide to quickly understanding the key parameters that need to be input to use the software. We have followed the original Python scripts as closely as possible, so it should be straightforward for users familiar with the origin Python scripts to migrate to the cpp code

### Documentation for individual executables in the suite

All the executables provide help and an example string that is accessible usint the -h or -? flag

### umisplit

#### Description

Takes pairs of R1 and R2 fastq files and generates a new fastq file with the UMI from the R2 file merged to the header line of the R1 file. It also separates the reads based on the plate well barcode part of the UMI. Each thread handles a different pair of input files.

#### Parameters

```
w96/umisplit parameters are: umisplit h?vft:m:N:o:b:l:q:
-h -? (display this message)
-v (Verbose mode)
-f (filter and discard reads with ambiguous UMI and barcode (default is to keep))
-z Compress the output as gzip files
-s The maximum size of the split file in KB

-l <Length of UMI (16):
-t <number of threads(1)>
-q <minimum quality of UMI base pair before changed to N (10)>-b <barcode file>
-m <maximum number of mismatches tolerated in barcode (0)>
-c <the file to output the umicounts to - default is UMIcounts.bin in the output directory>
-N <maximum number of ambiguous base pairs tolerated in barcode (0)>
-o <Output Directory>
<R1file.1> <R2file.1>..<R1file.N> <R2file.N>

Example:
umisplit -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz
```

### umimerge_parallel: 

#### Description
Takes sam files organized by plate wells and counts the reads mapped to each gene. The original gene level filtering scheme used by the Python script (i.e. reads with the same barcode and same gene mapping are filtered) is supported by the -g option. Otherwise position level filtering is used where reads with the same barcode that map to the same position in the gene are filtered which is the scheme from the original UMI paper. Each thread handles a different sam file.

#### Parameters

```
w96/umimerge_parallel parameters are: umimerge_parallel v?hfn:go:t:i:s:e:m:c:b:a:p:
-h -?  (display this message)
-v (Verbose mode)
-f Merge filtered SAM files (*.saf) instead of full SAM files
-s <sym2ref file>
-e <ercc_fasta file>
-b <barcode_file>
-a <aligns directory>
-o <dge_dir(counts directory)>
-c <binary file with umicounts - default is UMIcounts.bin in the input SAM directory>
-t <number of threads (1)>
-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... 

Required params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir

Example:

umimerge_parallel -i RNAseq_20150409 -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -a Aligns -o Counts -t 4 -p 0

```
### umimerge_filter:

umimerge_filter is experimental. It takes samfile output for an aligned read and writes a 16 bit hashvalue for the gene and map position.

```
umimerge_filter parameters are: umimerge_filter vh?i:g:n:e:b:o:
-h -?  (display this message)
-v (Verbose mode)
-g Filter identical UMIs that map to same gene
-i <sample_id>
-s <sym2ref file>
-e <ercc_fasta file>
-b <barcode_file>
-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... 

umimerge_filter is experimental. It takes samfile output for an aligned read and writes a 16 bit hashvalue for the gene and map position. Required params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir

Example:

cat myFile.sam | umimerge_filter  -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -p 0 > myOutputFile.saf
 
```

## Documentation of scripts provided

In addition to these scripts, the testScripts directory contains all the scripts used in the benchmarking in the manuscript that describes this work

### build.sh

A build script to construct the Docker container locally

### multibwa.sh 

A bash script to search a directory for sam files and run bwa in parallel. The input are files split into directories corresponding to the plate well barcode. 

###  start.sh

Used in the Docker container 

### fast_run-alignment-analysis.sh

See section on customizing shell scripts.

### docker_fast_run-alignment-analysis.sh 

See section on customizing shell scripts.

### runExample.sh

Downloads example files and runs example locally. The user must first compile the code before running.

### runDockerExample.sh
Downloads example files and runs the example using Docker. The user must install Docker and start Docker before running the script.

### download.sh
A bash script to download files - includes support for files on google drives.

## Summary of implementation changes
The major changes from the Python scripts are:

1. Use of unique numerical translation of UMI sequence and mapping position for hash to determine whether the read is a duplicate. This allows for deduplication of reads that map to the same or similar position rather than reads that map to the same gene.

2. Splitting reads by wells to increase speed, decrease memory requirements and allow for parallel processing of sam files. 

3. Use of a lookup table to detemine whether an ambiguous sequence maps uniquely to a known barcode

4. Slightly faster bwa alignment using UNIX threads to parallelize at the application level rather than relying on bwa's built-in parallelization.




