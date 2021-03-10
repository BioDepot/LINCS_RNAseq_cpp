#include <string>      
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <glob.h>

//default definitions - can be overridden by passing variables to make

#ifndef MAX_EDIT_DISTANCE
 #define MAX_EDIT_DISTANCE 1
#endif
#ifndef MAX_BEST
 #define MAX_BEST 20
#endif
#ifndef NWELLS
 #define NWELLS 96
#endif
#ifndef GLOB_TILDE
 #define GLOB_TILDE 0
#endif
#define  UMISIZE 10
#define  NUMIS 1048576 //4**10 - for 10 base UMI
#define SAMLINESIZE 1024
#define MAXLINESIZE 1024
//using namespace std;
uint32_t hashCode4(std::string &sequence);
bool ambigCheck(std::string &sequence,uint32_t &code);
std::string decodeId(uint32_t id, int size);
bool polyACheck(std::string &sequence);

void readRefseq(std::string filename, std::unordered_map<std::string,std::string> &refseqToGene,std::vector<std::string> &geneList);
std::string  readWells(std::string filename, std::unordered_map<std::string,uint32_t> &wellToIndex, std::vector<std::string> &wellList, std::vector<std::string> &shortWellIds);
std::string  readWells(std::string filename, std::unordered_map<std::string,uint32_t> &wellToIndex, std::vector<std::string> &wellList);
void readERCC(std::string filename, std::vector<std::string> &erccList);
void readAlignedFiles(std::string aligned_dir, std::vector<std::string> &alignedFiles, std::string suffix, std::string well, bool mixtureOfWells);
void splitStr(char *cstr,const char *delim, std::vector<std::string> &items);
void splitStr(std::string str,const char *delim, std::vector<std::string> &items);
std::string splitStrIndex(std::string str,const char *delim, int index);
bool multiGeneHit(std::vector<std::string> &best_list, std::string gene, const std::unordered_map<std::string,std::string> &refseqToGene);
void bitPrint(FILE *fp, uint64_t value,uint8_t valuebits);
uint8_t minBitsToRepresent(uint32_t number);
uint8_t bitsNeeded(uint32_t ncategories, uint8_t umiSize, uint8_t nbinbits, bool markMultiHits);
bool checkMultiHit(std::string alignedId, uint32_t pos, uint8_t nbinbits, uint8_t binsizebits, uint32_t posMask, std::string matchString, std::unordered_map<std::string,std::string> &refseqToGene, bool sameGeneHitNotMultiHit);

template <class T> T encodeMapping(uint32_t category, uint32_t umiIndex, uint32_t pos,uint8_t nbinbits, uint8_t binsizebits, uint8_t umibits, uint32_t posMask, T leftBitMask, bool multiHit, bool markMultiHits);
//no umi version 

template <class T> bool samToCategory(T &category, uint32_t &umiIndex, uint32_t &pos, bool &multiHit, std::string &alignedId, char *fullLine, uint8_t barcodeSize,uint8_t umiSize, std::unordered_map<std::string,std::string> &refseqToGene, std::unordered_map<std::string,uint32_t> &erccToIndex, std::unordered_map<std::string,uint32_t> &geneToIndex, uint32_t posMask, uint8_t binsizebits, uint8_t nbinbits, bool multiHits, bool sameGeneHitNotMultiHit, bool properPairs, bool &unassigned, bool &nonRefseq);

template <class T>class umipanel{
	//T is used to choose between unsigned char for 96 weill and uint16_t for 384 wells
	public:
	 T *hash;  //for barcode
	 std::vector<std::string> sequences; //contains panel of barcodes 
	 std::vector<std::string> wells;     //contains well of barcodes 
 	std::string panelID;
	 uint32_t nBarcodes=0;
	 uint32_t barcodeSize=0;
	 uint32_t hashSize=0;
	 uint32_t mismatchTol; //how many basepairs difference for matching panel
	 uint32_t NTol; //how many basepairs in query sequence
	 
	 umipanel(){
			nBarcodes=0;
			barcodeSize=0;
		}	
 	umipanel (std::string fileName,int _mismatchTol,int _NTol){
			//read in panel
			mismatchTol=_mismatchTol;
			NTol=_NTol;
 		std::string line,name,well,sequence;
   std::ifstream inFile(fileName,std::ifstream::in);
   if(!inFile){
				fprintf(stderr,"unable to read in barcode file %s\n",fileName.c_str());
				exit(EXIT_FAILURE);
			}
   while(getline(inFile,line)){
    std::istringstream iss(line);
    iss >> panelID; iss >> well;iss >> sequence;
    wells.push_back(well);
    sequences.push_back(sequence);
    if(!barcodeSize) barcodeSize=sequences[0].size();
    nBarcodes++;
 	 }
   inFile.close();
   hashSize=5;
   for(int i=1;i<barcodeSize;i++){
				hashSize*=5;
			}
   hash=new T[hashSize];
			memset(hash,0,hashSize*sizeof(T));
   fillHash();
	 }
	 ~umipanel(){
			if(hash)delete[] hash;
		}	
	 void fillHash(){
			char bp[5]={'N','A','C','G','T'};
			char *seq=new char[barcodeSize];
			char *indices=new char[barcodeSize];
			memset(indices,0,barcodeSize);
			hash[0]=0;
			for(int i=0;i<barcodeSize;i++)
			 seq[i]='N';
   for(int i=1;i<hashSize;i++){
				int numberofNs=0;
				int divisor=5;
				seq[0]=bp[i%divisor];
				int j=1;
				while(i%divisor == 0 && j<barcodeSize){
					indices[j]=(indices[j]+1)%5;
					seq[j]=bp[indices[j]];
					if(!seq[j])numberofNs++;
				 divisor*=5;
				 j++;
				}
				hash[i]=bestMatch(sequences,seq,nBarcodes,mismatchTol,NTol)+1;		
			}
			delete[] seq;	
			delete[] indices;	
		}
		uint32_t bestMatch(const char *query){
			return (uint32_t) hash[hashCode(query)];
		}		 
		uint32_t hashCode(const char *sequence){
			uint32_t code =0;
			int k=1;
			for (int i=0;i<barcodeSize;i++){
				switch (sequence[i]){
					case 'A':
					 code+=k;
					break;
					case 'C':
					 code+=2*k;
					 break;
					case 'G':
					 code+=3*k;
					 break;
				 case 'T':
				  code+=4*k;
				  break;  
				}
				k*=5;				
			}
			return code;	
		}
	uint32_t bestMatch (std::vector<std::string> &panelSeqs, char *query,int nPanelSeqs,int mismatchTol,int NTol){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		//we will initialize it anyway to get a meaningfull absolute distance
		int bestIndex=0,nBest=1;
		int maxDist=0;
		int nNs=0;
  for (uint32_t i = 0; i < barcodeSize; i++){
	 	if(query[i] != 'N' && panelSeqs[0][i] !=query[i]){
    maxDist++;
	 	}
	 	else if (query[i] == 'N') nNs++;
	 	if(nNs > NTol)return -1;
	 }
		for(int i=1;i<nPanelSeqs ;i++){
   int dist=0;
   for (uint32_t j = 0; j < barcodeSize && dist<=maxDist; j++){
				if(query[j] != 'N' && panelSeqs[i][j] !=query[j]){
     dist++;
				}	
			}	
		 if(maxDist == dist){				
				//if(!dist)return -1; //dupe found
				nBest++;
			}
			else if(dist < maxDist){
				nBest=1;
				bestIndex=i;
				maxDist=dist;
			}			
		}
		if(maxDist <=mismatchTol && nBest==1){
			return bestIndex;
		}
		return -1;	
	}			
};	
class MapPosition{
	public:
	std::vector <std::string> gene; 
	std::vector <int> position;
	MapPosition(){
		clear();
	}	
	bool insert(std::string _gene,int _position){
		for(int i=0;i<gene.size();i++){
	  if(_gene==gene[i] && _position==position[i]){
				return 0;
			}
		}
	 gene.push_back(_gene);
	 position.push_back(_position);
	 return 1;	
	}
	void clear(){
	 gene.clear();
	 position.clear();
 }
 bool multiGeneHit(std::string _gene, std::unordered_map<std::string,std::string> &refseqToGene){
	 for(int i=1;i<gene.size();i++){
	  if(!refseqToGene.count(gene[i]) || _gene != refseqToGene[gene[i]])return 1;
		}
	 return 0; 
 }	 	
};
uint32_t hashCode4(std::string &sequence){
		uint32_t code =0;
		int k=1;
		for (int i=0;i<sequence.size();i++){
		switch (sequence[i]){
				case 'C':
				 code+=k;
					 break;
					case 'G':
					 code+=2*k;
					 break;
				 case 'T':
				  code+=3*k;
				 break;  
				}
				k*=4;				
			}
	return code;	
}	
bool ambigCheck(std::string &sequence,uint32_t &code){	
	int k=1;
	code=0;
 for(int i=0;i<sequence.size();i++){
		switch (sequence[i]){
				 case 'N':
				  return 1;	
				 case 'C':
				  code+=k;
					break;
					case 'G':
					 code+=2*k;
					break;
				 case 'T':
				  code+=3*k;
				 break;  
				}
				k*=4;				
	}
	return 0;
}
std::string decodeId(uint32_t id, int size){
 std::string sequence;
 char bp[4]={'A','C','G','T'};
 for(int i=0;i<size;i++){
	 uint32_t index= (id >> 2*i) & 0x03;
  sequence.push_back(bp[index]);
	}
	return sequence;
	
}
bool multiGeneHit(std::vector<std::string> &best_list, std::string gene, std::unordered_map<std::string,std::string> &refseqToGene){
	//want to check that at least on of the alternative genes in the list is not the same as the top assignment 
	for (int i=0;i<best_list.size();i++)
		if(!refseqToGene.count(best_list[i]) || gene != refseqToGene[best_list[i]])return 1;
	return 0; 
}

void splitStr(char *cstr,const char *delim, std::vector<std::string> &items){
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	items.resize(0);		
	while(p){
		items.push_back(std::string(p));
		p=strtok_r(0,delim,&save);
	}	
}
void splitStr(std::string str,const char *delim, std::vector<std::string> &items){
	//have to make a copy of str in this case
	if(str.size()<1024){
		char cstr[1024];
		strcpy(cstr,str.c_str());
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	 items.resize(0);
	 while(p){
	 	items.push_back(std::string(p));
		 p=strtok_r(0,delim,&save);
	 }
	}
	else{
		char *cstr=(char*) malloc(str.size()+1);
		strcpy(cstr,str.c_str());
		char *save;
	 char *p=strtok_r(cstr,delim,&save);
	 items.resize(0);
	 while(p){
	 	items.push_back(std::string(p));
	 	p=strtok_r(0,delim,&save);
	 }
	 free(cstr);	
	}	 	
}
std::string splitStrIndex(std::string str,const char *delim, int index){
	std::vector <std::string> items;
	splitStr(str,delim,items);
	if (!items.size()) return "";
	if (index < 0 ){
		index=items.size()+index;
	}	
	if (index >= 0 && index < items.size()) return items[index];
 return "";
}

std::string readWells(std::string filename, std::unordered_map<std::string,uint32_t> &wellToIndex, std::vector<std::string> &wellList, std::vector<std::string> &shortWellIds){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64],id[64],well[64],seq[64];
	uint32_t k=0;
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s %s",id,well,seq);
	 wellList.push_back(std::string(id)+"_"+std::string(well));
	 shortWellIds.push_back(well);
	 wellToIndex[std::string(well)]=k;
	 wellToIndex[std::string(id)+"_"+std::string(well)]=k++;
	}
	fclose(fp);
	return std::string(id);	
}
std::string readWells(std::string filename, std::unordered_map<std::string,uint32_t> &wellToIndex,std::vector<std::string> &wellList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64],id[64],well[64],seq[64];
	uint32_t k=0;
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s %s",id,well,seq);
	 wellList.push_back(std::string(id)+"_"+std::string(well));
	 wellToIndex[std::string(well)]=k;
	 wellToIndex[std::string(id)+"_"+std::string(well)]=k++;
	}
	fclose(fp);
	return std::string(id);	
}	
void readERCC(std::string filename, std::vector<std::string> &erccList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64]; //max line width is 50
	while(fgets(line,sizeof(line),fp)){
	 if(line[0] == '>'){ //fasta file - just want comment line without carat
			char temp[64];
			sscanf(line+1,"%s",temp);
	  erccList.push_back(std::string(temp));
		}
	}
	fclose(fp);	
}

void readRefseq(std::string filename, std::unordered_map<std::string,std::string> &refseqToGene, std::vector<std::string> &geneList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[1024]; //max line width is 843
	char gene[256],refseq[1024];
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s",gene,refseq);
	 geneList.push_back(std::string(gene));
	 //replace first comma with 0
	 char *save;
	 char *p=strtok_r(refseq,",",&save);
	 while(p){
	  refseqToGene[std::string(p)]=std::string(gene);
	  p=strtok_r(0,",",&save);
		}
	}
	fclose(fp);	
}	

bool polyACheck(std::string &sequence){
	if(sequence.size() < 20) return 0;
	const char *c=sequence.c_str()+sequence.size()-20;
	while(*c){
		if(*c != 'A') return 0;
		c++;
	}
	return 1;	
}

bool readCountsFile(const char *fileName, uint32_t *counts){
	FILE *fp =fopen(fileName,"r");
	fprintf(stderr,"opening %s\n",fileName);
	if(!fp){
		fprintf(stderr,"error opening %s\n",fileName);
		return(0);
	}
	if(!fread(counts,sizeof(uint32_t)*NUMIS,1,fp)){
		fprintf(stderr,"error reading %s\n",fileName);
		fclose(fp);
		return 0;
	}
	fprintf(stderr,"read %s\n",fileName);	
	fclose(fp);
	return(1);
}			
bool filterSAMoutput(){
	char buffer[SAMLINESIZE];
	memset(buffer,0,sizeof(buffer));
	while(fgets(buffer, SAMLINESIZE, stdin)){
		if (buffer[0] != '@'){
			fputs(buffer,stdout);
		} 
 }
 return 1;
}
uint8_t minBitsToRepresent(uint32_t number){
    uint8_t nbits=0;
    while (number > 0){
	 nbits++;
	 number =  number >> 1;	
	}
	return nbits;		
}
uint8_t bitsNeeded(uint32_t ncategories, uint8_t umiSize, uint8_t nbinbits, bool markMultiHits){
    return minBitsToRepresent(ncategories) + nbinbits + markMultiHits + 2*umiSize;
}
void bitPrint(FILE *fp, uint64_t value, uint8_t valuebits){
	const uint64_t one=1;
	for (int i=valuebits-1;i>=0;i--){
			fprintf(fp,"%d",(value >>i) & one);
			if(!(i % 4))fprintf(fp," ");
	}
	fprintf(fp,"\n");
	
}
bool checkMultiHit(std::string alignedId, uint32_t pos, uint8_t nbinbits, uint8_t binsizebits, uint32_t posMask, std::string matchString, std::unordered_map<std::string,std::string> &refseqToGene, bool sameGeneHitNotMultiHit){
	std::vector<std::string> split_loc;
	std::string best_hits_loc = splitStrIndex(matchString,":",-1);
	splitStr(best_hits_loc,";",split_loc);
	//with ERCC and chrM we check the alignedId directly, otherwise we first transform to refSeq in case we have two names mapping to same id
	if (alignedId.substr(0,4) == "ERCC" || alignedId.substr(0,4) == "chrM" || !refseqToGene.count(alignedId)){
		for (auto iter = split_loc.begin(); iter != split_loc.end(); ++iter){
			std::string geneString=splitStrIndex(*iter,",",0);
			if (geneString != alignedId) return 1;
		}
	}
	else{
		std::string gene = refseqToGene[alignedId];
		for (auto iter = split_loc.begin(); iter != split_loc.end(); ++iter){
			std::string geneString=splitStrIndex(*iter,",",0);
			if(!refseqToGene.count(geneString) || gene != refseqToGene[geneString]) return 1;
		}
	}
	if (sameGeneHitNotMultiHit) return 0;
	//if we made it through to here then we check that the positions are the same
	if (nbinbits > 0 && nbinbits > binsizebits){
		uint32_t bin = pos >> binsizebits & posMask;
		for (auto iter = split_loc.begin(); iter != split_loc.end(); ++iter){
			std::string posString=splitStrIndex(*iter,",",1);
			if (posString=="") return 1;
			try {
				int splitPos = std::stoi(posString);
				if (splitPos == 0) return 1;
				if (splitPos > 0){
					uint32_t splitBin = (uint) splitPos >> binsizebits & posMask;
					if (splitBin != bin) return 1;
				}				
			}
			catch (...){
			    return 1;	
			}

		}		
	}
	return 0;
}


template <class T> T encodeMapping(uint32_t category, uint32_t umiIndex, uint32_t pos,uint8_t nbinbits, uint8_t binsizebits, uint8_t umibits, uint32_t posMask, T leftBitMask, bool multiHit, bool markMultiHits){
	T code=category;
//	bitPrint(stderr,code,64);
	if (umibits) code=(code << umibits) | umiIndex;
//	bitPrint(stderr,umiIndex,64);
//	bitPrint(stderr,code,64);
	if(nbinbits){
//		fprintf(stderr,"%d\n",pos);
		T positionCode=(pos >> binsizebits) & posMask;
//		bitPrint(stderr,positionCode,64);
		code = (code << nbinbits) | positionCode;
	}
//	bitPrint(stderr,code,64);	
	if (multiHit && markMultiHits) code = code | leftBitMask;
//	bitPrint(stderr,code,64);
	return code;
}

template <class T> bool samToCategory(T &category, uint32_t &umiIndex, uint32_t &pos, bool &multiHit, std::string &alignedId, char *fullLine, uint8_t barcodeSize,uint8_t umiSize, std::unordered_map<std::string,std::string> &refseqToGene, std::unordered_map<std::string,uint32_t> &erccToIndex, std::unordered_map<std::string,uint32_t> &geneToIndex, uint32_t posMask, uint8_t binsizebits, uint8_t nbinbits, bool multiHits, bool sameGeneHitNotMultiHit, bool properPairs, bool &unassigned, bool &nonRefseq){
	if(fullLine[0] == '@') return 0;
	const uint32_t geneListSize = geneToIndex.size();
	const uint32_t erccListSize = erccToIndex.size();
	//add 1 for chrM
	const uint32_t ncategories = geneListSize + erccListSize +1;
	unassigned=1;
	nonRefseq=0;
   
	//remove \n if it exists
	fullLine[strcspn(fullLine, "\r\n")] = 0;
	T code=0;
	std::vector <std::string> items;
	std::vector <std::string> tempItems;
	bool printFlag=0;
	splitStr(fullLine," \t",items);
	splitStr(items[0],":",tempItems);
	std::string fullBarcode,barcode;
	//umis then check for ambiguity in barcode
	if (barcodeSize || umiSize){
		//then we check for barcodes which are encoded here
		std::string fullBarcode=tempItems[tempItems.size()-1];
		if (umiSize){
			barcode=fullBarcode.substr(barcodeSize,umiSize);
		    //skip if ambiguous barcode and get unique barcode index from sequence
			if(ambigCheck(barcode,umiIndex)) return 0;
		}
	}
	//if it is not barcoded/UMI and properPairs is set that means we are looking at paired ends and will only count matched reads
	else if (properPairs){
		uint32_t flag=(uint32_t)stoi(items[0]);
		uint32_t mask=2;
		if (!(flag & mask)) return 0;
	}
	alignedId=std::string(items[2]);
	pos = stoi(items[3]);
	if(alignedId == "*") return 0;
	std::string read=items[9];
	//kallisto SAM does not have fields 10 if not aligned, does not have field 12 and 13 if extra
	int edit_dist=0,best_hits=0;
	
	if(items.size() == 14){
		best_hits=stoi(splitStrIndex(items[12],":",-1));
		multiHit = checkMultiHit(alignedId, pos, nbinbits, binsizebits, posMask, items[13],refseqToGene,sameGeneHitNotMultiHit);
	   //kallisto multihit
	}
	else if (items.size() > 19) {
		if  (items.size() > 12 && !items[12].compare(0,2,"X0")) edit_dist=stoi(splitStrIndex(items[12],":",-1));
		if  (items.size() > 13 && !items[13].compare(0,2,"X1")) best_hits=stoi(splitStrIndex(items[13],":",-1));
		//the original script does skip this read if any of these are true
		if(edit_dist > MAX_EDIT_DISTANCE || best_hits > MAX_BEST || polyACheck(read)) return 0;		
		if (best_hits >= 1 && items.size() > 19){
			multiHit = checkMultiHit(alignedId, pos, nbinbits, binsizebits, posMask, items[19],refseqToGene,sameGeneHitNotMultiHit);
		}
	}
	if (alignedId.substr(0,4) == "ERCC"){
		 category= geneListSize + erccToIndex[alignedId];
	}
	else if (alignedId.substr(0,4) == "chrM"){
		 category = ncategories-1;
	}
	else if (refseqToGene.count(alignedId)) category = geneToIndex[refseqToGene[alignedId]];
	else{
		nonRefseq=1;
		category=0;
	}
	return 1;
}	


