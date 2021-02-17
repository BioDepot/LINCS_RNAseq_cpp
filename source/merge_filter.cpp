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

extern "C" {
 #include "optparse.h"  
}
 //each code is associated with at least gene name or ercc spikein ID or chrM (mitochondrial) - the sum of possible categories is the nCategories
 //The unique encoding at the gene levels is barcodeIndex*nCategories+categoryID (geneID number if in refSeq or geneID+ERCC id if ERCC or last index if chrM
 //To add  (optional positional information we need to know the max number of bin and the bin size - for convenience we define this using number of bits
 //if only using coding genes or interested in unique hits (20203 coding + 17871 non-coding) then can use 16 bits
 //Otherwise need another bit for multigene or 17 bits
 //For titin 300K coding region or 19 bits -> 36 bits minimum for transcripts
 //barcodes are usually 6 bit -> 42 bit minimum


template <class T> void merge_filter(vector<string> &erccList, vector<string> &geneList, unordered_map<string,string> refseq_to_gene,unordered_map<string,unsigned int>ercc_to_index,unordered_map<string,unsigned int>gene_to_index, string inputFile, string outputFile){
    const T nCategories=geneList.size()+erccList.size()+1; //add1 for chrM
    const T leftBit = 1 << (sizeof(T)-1);
    
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
	    if(fullLine[0] == '@') continue;
		//remove \n if it exists
		T code;
	    fullLine[strcspn(fullLine, "\r\n")] = 0;
	    vector <string> items;
	    vector <string> tempItems;
	    bool printFlag=0;
	    splitStr(fullLine," \t",items);
	    splitStr(items[0],":",tempItems);
	    string aligned_id=items[2];
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
			code=erccIndex+geneList.size();
			if(best_hits_list.size()) {
				code=code | leftBit;	
			}
		}
		else if (aligned_id.substr(0,4) == "chrM"){
			code=nCategories-1;
			if(best_hits_list.size()){
				code=code | leftBit;	
			}
		}
		else if (refseq_to_gene.count(aligned_id)){
			string gene = refseq_to_gene[aligned_id];
			const unsigned int geneIndex =gene_to_index[gene];
			code=geneIndex;
			bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
			if(multiGene){
				code=code | leftBit;	
			}
		}
		else{
			//unknown
			continue;
		}
		fwrite(&code,1,sizeof(code),ofp);
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
 
 string errmsg="umimerge_filter vh?i:g:n:e:b:o:\n-h -?  (display this message)\n-v (Verbose mode)\n-g Filter identical UMIs that map to same gene\n-i <sample_id>\n-s <sym2ref file>\n-e <ercc_fasta file>\n Required params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir\n\nExample:\n\ncat myFile.sam | umimerge_filter  -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -p 0 > myOutputFile.saf\n";
 while ((opt = optparse(&options, "v?hgi:o:n:s:e:b:p:")) != -1) {
  switch (opt){
			case 'v':
			 verbose=1;
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
 

 uint8_t nbits=bitsNeeded(geneList.size(),erccList.size(),umiLength,binSize,keepMultiHits);
 if (nbits <= 16) merge_filter <uint16_t> (erccList,geneList,refseq_to_gene,ercc_to_index,gene_to_index,inputFile,outputFile);
 else if (nbits <= 32 ) merge_filter <uint32_t> (erccList,geneList,refseq_to_gene,ercc_to_index,gene_to_index,inputFile,outputFile);
 else merge_filter <uint64_t> (erccList,geneList,refseq_to_gene,ercc_to_index,gene_to_index,inputFile,outputFile);
 
 return 1;		
	 
}
