#include <zlib.h>  
#include <stdio.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <unordered_map> 
#include <unordered_set> 
#include <set> 
#include <vector> 
#include "kseq.h"
#include <glob.h>
#include "umitools.hpp"

//presently hardcoded for UMI size 16 and well barcode size 6
//to generalize have to generalize the hash function as well

 //barcodes are mapped as follows
 //for each umi convert obtain barcodeIndex - 2 bit conversion
 //00 A 01 C 10 G 11 T
 //N is present do not count - end search
 //ambiguity if it is to be handled is handled at the split/demux stage
 //This
 
 //each umi is associated with at least gene name or ercc spikein ID or chrM (mitochondrial) - the sum of possible categories is the nCategories
 //The unique encoding at the gene levels is barcodeIndex*nCategories+categoryID (geneID number if in refseq or geneID+ERCC id if ERCC or last index if chrM
 //This umi gives gene level unique encoding
 //To add positional information we need to know the max number of bin and the bin size - for convenience we define this using number of bits
 //First we calculate the bin by right shifting bin size bits and filtering through a mask of maxBins bits
 //left shift the unique encoding of the gene levels by maxBin bits and or it with the encoded bin 
 //The results should look like MxxxCCCCCCCCCCCUUUUUUUUUUUPPPPPPPPP   where the C's are the bits encoding the category, P are the bits encoding the position, U +bits for UMI and M the leftmost bit encoding whether there is a multihit
 //leftmost encoding avoids multiplication and modulo to get categories which unlike U are not strict powers of 2. P is not either but we leave so much room for the long genes that we have to assign a limit to the largest size whereas the number of genes/categories is easily counted 
 

extern "C" {
 #include "optparse.h"  
}

template <class T> void merge_filter(std::vector<std::string> &erccList ,std::vector<std::string> &geneList, std::unordered_map<std::string,std::string> refseqToGene,std::unordered_map<std::string,unsigned int>erccToIndex,std::unordered_map<std::string,unsigned int>geneToIndex, uint8_t binsizebits, uint8_t nbinbits, std::string inputFile, std::string outputFile, uint8_t barcodeSize, uint8_t umiSize, bool properPairs, bool markMultiHits, bool sameGeneHitNotMultiHit){
	
	const uint32_t geneListSize = geneList.size();
	const uint32_t erccListSize = erccList.size();
	const  uint32_t ncategories=geneListSize +erccListSize+1; //add1 for chrM - needs to be T type because of left bit shifts
	const T leftBitMask = 1ULL << (sizeof(T)-1);
	
	//eg nbinbits = 2
	//100 shift 2
	//11 subtract 1
	const uint32_t posMask=(1 << nbinbits) -1;
	
	char fullLine[MAXLINESIZE];
	memset(fullLine,0,sizeof(fullLine));
	FILE *ifp,*ofp;
	if (inputFile != "") ifp=fopen(inputFile.c_str(),"r");
	else ifp=stdin;
	if (outputFile != "") ofp=fopen(outputFile.c_str(),"w");
	else ofp=stdout;
	while(fgets(fullLine, MAXLINESIZE, ifp)){
		if(fullLine[0] == '@') continue;
		T category=0;
		uint32_t umiIndex=0,pos;
		bool multiHit=0, assignedFlag=0, nonRefseqFlag=0;
		std::string alignedId;
		if (samToCategory(category,umiIndex,pos,multiHit,alignedId,fullLine,barcodeSize,umiSize,refseqToGene,erccToIndex,geneToIndex,posMask, binsizebits, nbinbits,markMultiHits, sameGeneHitNotMultiHit, properPairs,assignedFlag,nonRefseqFlag)){
			if (markMultiHits || !multiHit){
				T code=encodeMapping(category,umiIndex,pos, nbinbits,binsizebits, 2*umiSize, posMask,leftBitMask, multiHit, markMultiHits);
				fprintf(stderr,"%d\n",category);
				fwrite(&code,1,sizeof(code),ofp);
			}
		}
    }
	if (ifp && ifp != stdin) fclose(ifp);
	if (ofp && ofp != stdout) fclose(ofp);
}
int main(int argc, char *argv[]){
	bool markMultiHits=0,markNonRefseq=0,properPairs=0;
	uint8_t barcodeSize=6, umiSize=10, nbinbits=16, binsizebits=0;
	bool geneLevelFilter=0,filteredSAMfiles=0,mixtureOfWells=0,sameGeneHitNotMultiHit=0;
	std::string sampleId="",sym2ref="", ercc_fasta="", barcodes="", alignDir="", resultsDir="",countsFile="", inputFile="", outputFile="";
	int opt,verbose=0,nthreads=1;
	struct optparse options;
 optparse_init(&options, argv);	
	std::string errmsg="umimerge_parallel v?hPM:g:i:o:s:e:p:B:U: -?  (display this message)\n-v (Verbose mode)\n--s <sym2ref file>\n-e <ercc_fasta file>\n-i <inputFile> (default stdin)\n-o <outputFile> (default stdout)\n-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... \n"; 
	while ((opt = optparse(&options, "v?hPMSi:o:n:s:e:p:B:U:")) != -1) {
		switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'n':
			 nbinbits=atoi(options.optarg);
			 //set nbinbits - default is 16 bits otherwise  
			break;
			case 'i':
				inputFile=std::string(options.optarg);
			break; 
			case 'o':
				outputFile=std::string(options.optarg);
			break; 
			case 's':
				sym2ref=std::string(options.optarg);
			break;
			case 'e':
				ercc_fasta=std::string(options.optarg);
			break;   
			case 'p':		 
				binsizebits=atoi(options.optarg);
			break;
			case 'P':		 
				properPairs=1;
				//for paired ends 
			break;
			case  'M':
				markMultiHits=1;
			break;
			case  'B':
				barcodeSize=atoi(options.optarg);;
			break;
			case  'U':
				umiSize=atoi(options.optarg);
			break;
			case 'S':
				sameGeneHitNotMultiHit=1;
			break;
			case '?':
				fprintf(stderr, "%s parameters are: %s\n", argv[0], errmsg.c_str());
				exit(0);
			break;
			case 'h':
				fprintf(stderr, "%s parameters are: %s\n", argv[0], errmsg.c_str());
			exit(0);
			break; 
		}
	}
	if(sym2ref==""){
		fprintf(stderr,"Required parameter is-s sym2ref \n");
		exit(EXIT_FAILURE);
 }	
 
 std::unordered_map<std::string,std::string> refseqToGene;
 std::vector<std::string>erccList, geneList;
 std::unordered_map<std::string,uint32_t>erccToIndex;
 std::unordered_map<std::string,uint32_t>geneToIndex;
 std::vector<std::string> unknown_list;

 if (ercc_fasta != "") readERCC(ercc_fasta,erccList);
 if (sym2ref != "") readRefseq(sym2ref,refseqToGene ,geneList);
 for(int i=0;i<geneList.size();i++)
	geneToIndex[geneList[i]]=i;
 for(int i=0;i<erccList.size();i++)
	erccToIndex[erccList[i]]=i; 

uint8_t bitSize=bitsNeeded(erccList.size()+geneList.size()+1,umiSize, nbinbits, markMultiHits);
uint32_t posMask=(1 << nbinbits) -1;

 if (bitSize <= 16)  merge_filter<uint16_t>(erccList,geneList,refseqToGene,erccToIndex,geneToIndex,binsizebits,nbinbits,inputFile,outputFile,barcodeSize,umiSize,properPairs,markMultiHits,sameGeneHitNotMultiHit);
 else if (bitSize <= 32 )  merge_filter<uint32_t>(erccList,geneList,refseqToGene,erccToIndex,geneToIndex,binsizebits,nbinbits,inputFile,outputFile,barcodeSize,umiSize,properPairs,markMultiHits,sameGeneHitNotMultiHit);
 else merge_filter<uint64_t>(erccList,geneList,refseqToGene,erccToIndex,geneToIndex,binsizebits,nbinbits,inputFile,outputFile,barcodeSize,umiSize,properPairs,markMultiHits,sameGeneHitNotMultiHit);
 return 1;		
	 
}		
	 

