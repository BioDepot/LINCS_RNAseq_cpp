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
#define MAXLINESIZE 1024
//presently hardcoded for UMI size 16 and well barcode size 6
//to generalize have to generalize the hash function as well

 //barcodes are mapped as follows
 //for each umi convert obtain barcodeIndex - 2 bit conversion
 //00 A 01 C 10 G 11 T
 //N is present do not count - end search
 //ambiguity if it is to be handled is handled at the split/demux stage
 //This
 
 //each umi is associated with at least gene name or ercc spikein ID or chrM (mitochondrial) - the sum of possible categories is the nCategories
 //The unique encoding at the gene levels is barcodeIndex*nCategories+categoryID (geneID number if in refSeq or geneID+ERCC id if ERCC or last index if chrM
 //This umi gives gene level unique encoding
 //To add positional information we need to know the max number of bin and the bin size - for convenience we define this using number of bits
 //First we calculate the bin by right shifting bin size bits and filtering through a mask of maxBins bits
 //left shift the unique encoding of the gene levels by maxBin bits and or it with the encoded bin 
 //The results should look like CCCCCCCCCCCPPPPPPPPP   where the C's are the bits encoding the category/UMI and P are the bits encoding the position

 

extern "C" {
 #include "optparse.h"  
}
	
template <class T1,class T2> void merge_filter(vector<string> &erccList ,vector<string> &geneList, umipanel<T1,T2> &barcodePanel,  unordered_map<string,string> refseq_to_gene,unordered_map<string,unsigned int>ercc_to_index,unordered_map<string,unsigned int>gene_to_index, uint32_t posMask, int binSize, int nbins, bool geneLevelFilter, string inputFile, string outputFile){
		const uint64_t nCategories=geneList.size()+erccList.size()+1; //add1 for chrM
		const uint64_t bit63 = 1ULL << 63;
	 //The encoding here is slightly different we use the leftmost bit to indicate whether it is multigene or not
	 //This means that 1 fewer bit is available for encoding the hashes
  char fullLine[MAXLINESIZE];
  memset(fullLine,0,sizeof(fullLine));
  FILE *ifp,*ofp;
  if (inputFile != "") ifp=fopen(inputFile.c_str(),"r");
  else ifp=stdin;
  if (outputFile != "") ofp=fopen(outputFile.c_str(),"w");
  else ofp=stdout;
  while(fgets(fullLine, MAXLINESIZE, ifp)){
		 if(fullLine[0] != '@'){
				//remove \n if it exists
				fullLine[strcspn(fullLine, "\r\n")] = 0;
				uint64_t umicode=0;
				vector <string> items;
				vector <string> tempItems;
				bool printFlag=0;
				splitStr(fullLine," \t",items);
				splitStr(items[0],":",tempItems);
 	  string fullBarcode=tempItems[tempItems.size()-1];
	   string barcode=fullBarcode.substr(6,UMISIZE); //change if barcode size changes
	   
	   //ambigCheck
	   //skip if ambiguous barcode and get unique barcode index from sequence
	   unsigned int barcodeIndex=0;
	   if(ambigCheck(barcode,barcodeIndex))continue;
	   string aligned_id=items[2];
				unsigned pos = stoi(items[3]);
	   if(aligned_id == "*")continue;
	   string read=items[9];
	   int edit_dist=stoi(splitStrIndex(items[12],":",-1));
	   int best_hits=stoi(splitStrIndex(items[13],":",-1));
	   //the original script does skip this read if any of these are true
    if(edit_dist > MAX_EDIT_DISTANCE || best_hits > MAX_BEST || polyACheck(read)) continue;
    vector<string> best_hits_list;
    if(items.size() > 19){
					string best_hits_loc = splitStrIndex(items[19],":",-1);
	    vector<string> split_loc;
	    splitStr(best_hits_loc,";",split_loc);
					for (auto iter = split_loc.begin(); iter != split_loc.end(); ++iter){
					 string geneString=splitStrIndex(*iter,",",0);
					 if (geneString != "")
						 best_hits_list.push_back(geneString);			
					}
				}
				if (aligned_id.substr(0,4) == "ERCC"){
					const int erccIndex=ercc_to_index[aligned_id];
					umicode=barcodeIndex*nCategories+erccIndex+geneList.size();
					if(!geneLevelFilter){
						uint32_t positionCode=(pos >> binSize) & posMask; 
						umicode = (umicode << nbins) | positionCode;
					}
					if(best_hits_list.size()) {
						umicode=umicode | bit63;	
					}
				}
				else if 	(aligned_id.substr(0,4) == "chrM"){
					umicode=barcodeIndex*nCategories+nCategories-1;
					if(!geneLevelFilter){
						uint32_t positionCode=(pos >> binSize) & posMask; 
						umicode = (umicode << nbins) | positionCode;
					}
					if(best_hits_list.size()){
						umicode=umicode | bit63;	
					}
				}
				else if (refseq_to_gene.count(aligned_id)){
					string gene = refseq_to_gene[aligned_id];
					const unsigned int geneIndex =gene_to_index[gene];
					umicode=barcodeIndex*nCategories+geneIndex;
					if(!geneLevelFilter){
						uint32_t positionCode=(pos >> binSize) & posMask; 
						umicode = (umicode << nbins) | positionCode;
					}
					bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
					if(multiGene){
						umicode=umicode | bit63;	
					}
				}
				else{
					//unknown
					continue;
				}
				fwrite(&umicode,1,sizeof(umicode),ofp);
			}
		}
		if (ifp && ifp != stdin) fclose(ifp);
		if (ofp && ofp != stdout) fclose(ofp);
	}

using namespace std;   
int main(int argc, char *argv[]){
	int nbins=16;
	int binSize=0;
	uint32_t posMask=(1 << nbins+1) -1;
	bool geneLevelFilter=0;
	string sample_id="",sym2ref="", ercc_fasta="", barcodes="",inputFile="",outputFile="";
	int opt,verbose=0;
	struct optparse options;
 optparse_init(&options, argv);	
 
 string errmsg="umimerge_filter vh?i:g:n:e:b:o:\n-h -?  (display this message)\n-v (Verbose mode)\n-g Filter identical UMIs that map to same gene\n-i <sample_id>\n-s <sym2ref file>\n-e <ercc_fasta file>\n-b <barcode_file>\n-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... \n\numimerge_filter is experimental. It takes samfile output for an aligned read and writes a 16 bit hashvalue for the gene and map position. Required params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir\n\nExample:\n\ncat myFile.sam | umimerge_filter  -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -p 0 > myOutputFile.saf\n";
 while ((opt = optparse(&options, "v?hgi:o:n:s:e:b:p:")) != -1) {
  switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'n':
			 nbins=atoi(options.optarg);
			 posMask=(1 << nbins+1) -1;
			 //set nbins - default is 16 bits otherwise  
			break;
			case 'g':
			 //set gene level filtering - overrides other options
			 geneLevelFilter=1;
			break;
			case 'i':
    inputFile=string(options.optarg);
   break; 
   case 'o':
    outputFile=string(options.optarg);
   break;  
   case 's':
    sym2ref=string(options.optarg);
   break;
   case 'e':
    ercc_fasta=string(options.optarg);
   break; 
   case 'b':
    barcodes=string(options.optarg);
   break;
   case 'p':		 
    binSize=atoi(options.optarg);
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
 if(sym2ref=="" ||  ercc_fasta=="" || barcodes==""){
		fprintf(stderr,"Required params are -s sym2ref -e ercc_fasta -b barcodes\n");
		exit(EXIT_FAILURE);
	}	
 

 unordered_map<string,string> refseq_to_gene;
 vector<string>erccList, geneList;
 unordered_map<string,unsigned int>ercc_to_index;
 unordered_map<string,unsigned int>gene_to_index;
 vector<string> unknown_list;

 readERCC(ercc_fasta,erccList);
 readRefseq(sym2ref,refseq_to_gene ,geneList);
 
 for(int i=0;i<geneList.size();i++)
  gene_to_index[geneList[i]]=i;
 for(int i=0;i<erccList.size();i++)
  ercc_to_index[erccList[i]]=i;
 
 //outputfile names

 int mismatchTol=0;
 int NTol=0;
 
 #if NWELLS > 256
 umipanel<uint32_t,uint16_t> barcodePanel(barcodes.c_str(),mismatchTol,NTol);  		
	#else
	umipanel<uint32_t,unsigned char> barcodePanel(barcodes.c_str(),mismatchTol,NTol);  		
	#endif
	
 merge_filter(erccList,geneList,barcodePanel,refseq_to_gene,ercc_to_index,gene_to_index,posMask,binSize,nbins,geneLevelFilter,inputFile,outputFile);
 return 1;		
	 
}
