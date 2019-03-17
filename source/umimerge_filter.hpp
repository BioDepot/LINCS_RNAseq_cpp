#include <string>      
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>

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
using namespace std;
unsigned int hashCode4(string &sequence);
bool ambigCheck(string &sequence,unsigned int &code);
string decodeId(unsigned int id, int size);
bool polyACheck(string &sequence);

void readRefseq(string filename, unordered_map<string,string> &refseq_to_gene,vector<string> &geneList);
string readWells(string filename, unordered_map<string,unsigned int> &well_to_index,vector<string> &wellList);
void readERCC(string filename, vector<string> &erccList);
void splitStr(char *cstr,const char *delim, vector<string> &items);
void splitStr(string str,const char *delim, vector<string> &items);
string splitStrIndex(string str,const char *delim, int index);
bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene);
unsigned int sumCounts(unsigned int *counts,unsigned int size1);
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2);
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2);

template <class T1, class T2>class umipanel{
	//T1 is used to choose between 32 and 64 bit uint
	//T2 is used to choose between unsigned char for 96 weill and uint16_t for 384 wells
	public:
	 T2 *hash;  //for barcode
	 vector<string> sequences; //contains panel of barcodes 
	 vector<string> wells;     //contains well of barcodes 
 	string panelID;
	 unsigned int nBarcodes=0;
	 unsigned int barcodeSize=0;
	 unsigned int hashSize=0;
	 unsigned int mismatchTol; //how many basepairs difference for matching panel
	 unsigned int NTol; //how many basepairs in query sequence
	 
	 umipanel(){
			nBarcodes=0;
			barcodeSize=0;
		}	
 	umipanel (string fileName,int _mismatchTol,int _NTol){
			//read in panel
			mismatchTol=_mismatchTol;
			NTol=_NTol;
 		string line,name,well,sequence;
   ifstream inFile(fileName,ifstream::in);
   if(!inFile){
				fprintf(stderr,"unable to read in barcode file %s\n",fileName.c_str());
				exit(EXIT_FAILURE);
			}
   while(getline(inFile,line)){
    istringstream iss(line);
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
   hash=new T2[hashSize];
			memset(hash,0,hashSize*sizeof(T2));
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
		unsigned int bestMatch(const char *query){
			return (unsigned int) hash[hashCode(query)];
		}		 
		unsigned int hashCode(const char *sequence){
			unsigned int code =0;
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
	unsigned int bestMatch (vector<string> &panelSeqs, char *query,int nPanelSeqs,int mismatchTol,int NTol){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		//we will initialize it anyway to get a meaningfull absolute distance
		int bestIndex=0,nBest=1;
		int maxDist=0;
		int nNs=0;
  for (unsigned int i = 0; i < barcodeSize; i++){
	 	if(query[i] != 'N' && panelSeqs[0][i] !=query[i]){
    maxDist++;
	 	}
	 	else if (query[i] == 'N') nNs++;
	 	if(nNs > NTol)return -1;
	 }
		for(int i=1;i<nPanelSeqs ;i++){
   int dist=0;
   for (unsigned int j = 0; j < barcodeSize && dist<=maxDist; j++){
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
	vector <string> gene; 
	vector <int> position;
	MapPosition(){
		clear();
	}	
	bool insert(string _gene,int _position){
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
 bool multiGeneHit(string _gene, unordered_map<string,string> &refseq_to_gene){
	 for(int i=1;i<gene.size();i++){
	  if(!refseq_to_gene.count(gene[i]) || _gene != refseq_to_gene[gene[i]])return 1;
		}
	 return 0; 
 }	 	
};
unsigned int hashCode4(string &sequence){
		unsigned int code =0;
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
bool ambigCheck(string &sequence,unsigned int &code){	
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
string decodeId(unsigned int id, int size){
 string sequence;
 char bp[4]={'A','C','G','T'};
 for(int i=0;i<size;i++){
	 unsigned int index= (id >> 2*i) & 0x03;
  sequence.push_back(bp[index]);
	}
	return sequence;
	
}
bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene){
	//want to check that at least on of the alternative genes in the list is not the same as the top assignment 
	for(int i=0;i<best_list.size();i++)
	 if(!refseq_to_gene.count(best_list[i]) || gene != refseq_to_gene[best_list[i]])return 1;
	return 0; 
}	
void splitStr(char *cstr,const char *delim, vector<string> &items){
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	items.resize(0);		
	while(p){
		items.push_back(string(p));
		p=strtok_r(0,delim,&save);
	}	
}
void splitStr(string str,const char *delim, vector<string> &items){
	//have to make a copy of str in this case
	if(str.size()<1024){
		char cstr[1024];
		strcpy(cstr,str.c_str());
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	 items.resize(0);
	 while(p){
	 	items.push_back(string(p));
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
	 	items.push_back(string(p));
	 	p=strtok_r(0,delim,&save);
	 }
	 free(cstr);	
	}	 	
}
string splitStrIndex(string str,const char *delim, int index){
	vector <string> items;
	splitStr(str,delim,items);
	if (!items.size()) return "";
	if (index < 0 ){
		index=items.size()+index;
	}	
	if (index >= 0 && index < items.size()) return items[index];
 return "";
}

string readWells(string filename, unordered_map<string,unsigned int> &well_to_index,vector<string> &wellList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64],id[64],well[64],seq[64];
	unsigned int k=0;
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s %s",id,well,seq);
	 wellList.push_back(string(id)+"_"+string(well));
	 well_to_index[string(well)]=k;
	 well_to_index[string(id)+"_"+string(well)]=k++;
	}
	fclose(fp);
	return string(id);	
}
	
void readERCC(string filename, vector<string> &erccList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64]; //max line width is 50
	while(fgets(line,sizeof(line),fp)){
	 if(line[0] == '>'){ //fasta file - just want comment line without carat
			char temp[64];
			sscanf(line+1,"%s",temp);
	  erccList.push_back(string(temp));
		}
	}
	fclose(fp);	
}

void readRefseq(string filename, unordered_map<string,string> &refseq_to_gene, vector<string> &geneList){
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[1024]; //max line width is 843
	char gene[256],refseq[1024];
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s",gene,refseq);
	 geneList.push_back(string(gene));
	 //replace first comma with 0
	 char *save;
	 char *p=strtok_r(refseq,",",&save);
	 while(p){
	  refseq_to_gene[string(p)]=string(gene);
	  p=strtok_r(0,",",&save);
		}
	}
	fclose(fp);	
}	

bool polyACheck(string &sequence){
	if(sequence.size() < 20) return 0;
	const char *c=sequence.c_str()+sequence.size()-20;
	while(*c){
		if(*c != 'A') return 0;
		c++;
	}
	return 1;	
}
unsigned int sumCounts(unsigned int *counts,unsigned int size1){
	unsigned int sum=0;
	for(int i=0;i<size1;i++)
	  sum+=counts[i];
	return sum;
}	
		
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2){
	unsigned int sum=0;
	for(int i=0;i<size1;i++)
	 for(int j=0;j<size2;j++)
	  sum+=counts[i][j];
	return sum;
}
		
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2){
	unsigned int sum=0;
	for(int j=0;j<size2;j++)
	 sum+=counts[i][j];
	return sum;
}
bool readCountsFile(const char *fileName, unsigned int *counts){
	FILE *fp =fopen(fileName,"r");
	fprintf(stderr,"opening %s\n",fileName);
	if(!fp){
		fprintf(stderr,"error opening %s\n",fileName);
		return(0);
	}
	if(!fread(counts,sizeof(unsigned int)*NUMIS,1,fp)){
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
