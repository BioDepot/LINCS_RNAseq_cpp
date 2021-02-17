#include <zlib.h>  
#include <omp.h>
#include <stdio.h>
#include <string.h>
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
 
 //each umi is associated with at least gene name or ercc spikein ID or chrM (mitochondrial) - the sum of possible categories is the ncategories
 //The unique encoding at the gene levels is barcodeIndex*ncategories+categoryID (geneID number if in refseq or geneID+ERCC id if ERCC or last index if chrM
 //This umi gives gene level unique encoding
 //To add positional information we need to know the max number of bin and the bin size - for convenience we define this using number of bits
 //First we calculate the bin by right shifting bin size bits and filtering through a mask of maxBins bits
 //left shift the unique encoding of the gene levels by maxBin bits and or it with the encoded bin 
 //The results should look like CCCCCCCCCCCPPPPPPPPP   where the C's are the bits encoding the category/UMI and P are the bits encoding the position

 

extern "C" {
 #include "optparse/optparse.h"  
}

//some of the categories are holdovers from how the original counts were done	
template <class T> class Counts{
	private:
		//gene and well independent totals - but are dependent on number of threads
		//mm suffix indicates that multiHits with different geneposition are included - if nbinbits is zero bits then there is only one bin and all positions are treated the same within the gene for the purpose of multihit determination
		  
		uint32_t *total_reads, *assigned_reads, *aligned_reads,*assigned_unknown_reads, *assigned_unknown_umi, *total_reads_mm, *assigned_reads_mm, *aligned_reads_mm, *assigned_unknown_reads_mm, *assigned_unknown_umi_mm;
		
		//category dependent totals - also dependent on wells
		//each thread will work on a separate well so we don't need to partition for threads
		uint32_t **categoryTotal, **categoryUmi, **categoryTotal_mm, **categoryUmi_mm;

		//this hash is necessary to join the unknown hashes	 
		std::unordered_set<std::string>unknown_set;
		std::unordered_set<std::string> *unknown_set_w;
		std::unordered_set<std::string> unknown_set_mm;
		std::unordered_set<std::string> *unknown_set_w_mm;
		
		const uint32_t nthreads,erccListSize,geneListSize,ncategories,nwells;
		const uint8_t barcodeSize,umiSize;
		
		const bool markUnknowns,markMultiHits,mixtureOfWells,properReads;
		
		const std::vector<std::string> wellList;
		const std::vector<std::string> categoryList;
		const std::string alignsDir;
		const std::string sampleId;

		const std::unordered_map<std::string,std::string> refseqToGene;
		const std::unordered_map<std::string,uint32_t>erccToIndex;
		const std::unordered_map<std::string,uint32_t>geneToIndex;
		const uint32_t posMask;
		const uint8_t binsizebits, nbinbits, umibits;

		const T leftBitMask;
		
		uint32_t sumThreadCounts(uint32_t *counts){
			if(!counts) return 0;
			uint32_t sum=0;
			for(int i=0;i<nthreads;i++)
				sum+=counts[i];
			return sum;
		}	
		uint32_t sumCategoryCounts(uint32_t **counts, uint32_t offset,uint32_t categorySize){
			if(!counts || !categorySize ) return 0;
			uint32_t sum=0;
			for(int i=0;i<nwells;i++)
				for(int j=0;j<categorySize;j++)
					sum+=counts[i][j];
			return sum;
		}
		
		uint32_t sumWellCounts(uint32_t **counts,uint32_t wellIndex, uint32_t offset, uint32_t categorySize){
			if(!counts || !categorySize ) return 0;
			uint32_t sum=0;
			for(int j=offset;j<offset+categorySize;j++)
				sum+=counts[wellIndex][j];
			return sum;
		}
	
		void printWellSums(FILE *fp, std::string label,uint32_t **counts, uint32_t offset, uint32_t categorySize){		
			fprintf(fp,label.c_str()); 
			for(int i=0;i<nwells;i++){
				fprintf(fp,"\t%d",sumWellCounts(counts,i,offset,categorySize));
			}
		    fprintf(fp,"\n");
		}
			
		void readAlignedFiles(std::string alignDir, std::vector<std::string> &alignedFiles, std::string suffix, std::string well, bool 	mixtureOfWells){
			glob_t glob_result;
			std::string globString=alignDir+"/"+well+"/"+"*." + suffix;
			if (mixtureOfWells) globString=alignDir+"/*_"+well+"_"+"*."+suffix;
			glob(globString.c_str(),GLOB_TILDE,NULL,&glob_result);
			fprintf(stderr,"The glob string is %s\n",globString.c_str());
			for(int j=0; j<glob_result.gl_pathc; j++){
				alignedFiles.push_back(std::string(glob_result.gl_pathv[j]));
			}	
		}
		void printCounts (std::string filename,uint32_t categoryStart, uint32_t numprint, uint32_t **counts){
			if(!counts || !numprint) return;
			FILE *fp=0;
			if(!(fp=fopen(filename.c_str(),"w")))exit(EXIT_FAILURE);
			for(int i=0;i<nwells;i++){
				fprintf(fp,"\t%s",wellList[i].c_str());
			}
			fprintf(fp,"\n");
			for(int i=categoryStart;i<categoryStart+numprint;i++){
				fprintf(fp,"%s\t",categoryList[i].c_str());
				for(int j=0;j<nwells-1;j++){
					fprintf(fp,"%d\t",counts[j][i]);
				}
				fprintf(fp,"%d\n",counts[nwells-1][i]);
			}
			fclose(fp);
		}
		uint32_t decodeMappings (std::string inputFile){
			uint32_t nreads=0;
			FILE *fp=fopen(inputFile.c_str(),"r");
			T code;
			if(markMultiHits){
				while (fread(&code,sizeof(T),1,fp)){
					const T category= (code >> nbinbits);
					if (code >> (sizeof(T)-1)) categoryTotal_mm[category]++;
					else categoryTotal[category]++;	
					nreads++;
				}				
			}
			else{
				while (fread(&code,sizeof(T),1,fp)){
					const T category= (code >> nbinbits);
					categoryTotal[category]++;
					nreads++;	
				}
			}
			fclose(fp);
			return nreads;
		}
		uint32_t decodeMappingsUmi (std::string inputFile, uint8_t umipositionBits,std::unordered_set<T> &umi_seen, std::unordered_set<T> &umi_seen_mm){
			uint32_t nreads=0;
			T umicode;
			//read in each file
			FILE *fp=fopen(inputFile.c_str(),"r");
			if (markMultiHits){
				while (fread(&umicode,sizeof(T),1,fp)){
					const T category= (umicode >> umipositionBits);
					categoryTotal_mm[category]++;
					if ( umicode >> (sizeof(T)-1)){
						//multiHit
						if(!umi_seen_mm.count(umicode)){
							umi_seen_mm.insert(umicode);
							categoryUmi_mm[category]++;
						}
					}
					else{
						categoryTotal[category]++;
						if(!umi_seen.count(umicode)){
							umi_seen.insert(umicode);
							categoryUmi[category]++;
						}
					}
					nreads++;				
				}	
			}
			else{
				while (fread(&umicode,sizeof(T),1,fp)){
					const T category= (umicode >> umipositionBits);
					categoryTotal[category]++;
					if(!umi_seen.count(umicode)){
						umi_seen.insert(umicode);
						categoryUmi[category]++;	
					}
					nreads++;				
				}
			}
			fclose(fp);
			return nreads;
		}
	 	void updateCounts(T category, T code, bool multiHit, std::string *unknown, std::unordered_set<T>umi_seen, std::unordered_set<T>umi_seen_mm, uint32_t tid){
			if (multiHit){
				if (!markMultiHits) return;
				aligned_reads_mm[tid]++;
				if (unknown){
					if(markUnknowns){
						assigned_unknown_reads_mm[tid]++;
						if (!unknown_set_w_mm[tid].count(*unknown)){
							unknown_set_w_mm[tid].insert(*unknown);
						}
					}	
					return;
				}
				assigned_reads_mm[tid]++;
				categoryTotal_mm[tid][category]++;
				if(!umi_seen_mm.count(code)){
					umi_seen_mm.insert(code);
					categoryUmi_mm[tid][category]++;
				}		
			}	
			else{
				aligned_reads[tid]++;
				if (unknown){
					if(markUnknowns){
						assigned_unknown_reads[tid]++;
						if (!unknown_set_w[tid].count(*unknown)){
							unknown_set_w[tid].insert(*unknown);
						}
					}	
					return;
				}
				aligned_reads_mm[tid]++;
				categoryTotal[tid][category]++;
				if(!umi_seen.count(code)){
					umi_seen.insert(code);
					categoryUmi[tid][category]++;
				}	
			}				
		}
		void updateCounts(T category, bool multiHit,std::string *unknown, uint32_t tid){
			if (multiHit){
				if (!markMultiHits) return;
				aligned_reads_mm[tid]++;
				if (unknown){
					if(markUnknowns) assigned_unknown_reads_mm[tid]++;
					return;
				}
				assigned_reads_mm[tid]++;
				categoryTotal_mm[tid][category]++;	
			}	
			else{
				aligned_reads[tid]++;
				if (unknown){
					if(markUnknowns) assigned_unknown_reads[tid]++;
					return;
				}
				assigned_reads[tid]++;
				categoryTotal[tid][category]++;
			}				
		}
		void merge_sam(){		
			#pragma omp parallel for num_threads (nthreads) schedule (dynamic)
			for	(int wellIndex=0; wellIndex<nwells; wellIndex++){
				const uint32_t tid=omp_get_thread_num();
				const std::string well= wellList[wellIndex];
				std::vector<std::string> inputFiles;
				std::unordered_set<T> umi_seen, umi_seen_mm;		
				readAlignedFiles(alignsDir,inputFiles,"sam",well,mixtureOfWells);
				for(int i=0;i<inputFiles.size();i++){
					char fullLine[MAXLINESIZE];
					memset(fullLine,0,sizeof(fullLine));
					FILE *fp=fopen(inputFiles[i].c_str(),"r");
					while(fgets(fullLine, MAXLINESIZE, fp)){
						bool assignedFlag=0, unknownFlag=0;
						T category=0;
						uint32_t umiIndex=0,pos;
						bool multiHit=0;
						std::string alignedId;
						if (samToCategory<T>(category,umiIndex,pos,multiHit,alignedId,fullLine,barcodeSize,umiSize,refseqToGene,erccToIndex,geneToIndex,posMask, binsizebits, nbinbits,markMultiHits, properReads,assignedFlag,unknownFlag)){
							std::string *unknown=0;
							if (unknownFlag && markUnknowns) unknown=&alignedId;
							if (umiSize > 0){
								T code=encodeMapping(category,umiIndex,pos, nbinbits,binsizebits, umibits, posMask,leftBitMask, multiHit, markMultiHits);
								updateCounts(category,code,multiHit,unknown,umi_seen, umi_seen_mm,tid);
							}
							else{
								updateCounts(category,multiHit,unknown,tid);
							}
						}
					}
				}		
			}
		}	
	//the generalized version of the saf reader - no unknown_reads are kept
		void merge_saf(){	
			const T ncategories=geneListSize +erccListSize+1; //add1 for chrM
			const uint32_t nWells = wellList.size() ? wellList.size() : 1;
			#pragma omp parallel for num_threads (nthreads) schedule (dynamic)
			for	(int wellIndex=0; wellIndex<nwells; wellIndex++){
				const uint32_t tid=omp_get_thread_num();
				const std::string well= wellList.size()? wellList[wellIndex] : "all";
				std::vector<std::string> inputFiles;
				readAlignedFiles(alignsDir,inputFiles,"saf",well,mixtureOfWells);
				if (umibits){
					std::unordered_set<T>umi_seen;
					std::unordered_set<T>umi_seen_mm;
					for(int i=0;i<inputFiles.size();i++){
						aligned_reads[tid]+=decodeMappingsUmi(inputFiles[i].c_str(),nbinbits+umibits,umi_seen, umi_seen_mm);
	                }
	            }
	            else{
					for(int i=0;i<inputFiles.size();i++){					
						aligned_reads[tid]+=decodeMappings(inputFiles[i].c_str());
	                }
				}
			}
		}
	
	public:
		Counts(uint32_t in_nthreads, uint32_t maxErccListSize, uint32_t maxGeneListSize, bool in_markMultiHits, bool in_mixtureOfWells, bool in_markUnknowns,bool in_properReads,uint8_t in_barcodeSize, uint8_t in_umiSize, uint8_t in_nbinbits, uint8_t in_binsizebits, std::string in_alignsDir, std::string in_sampleId, std::vector<std::string>& in_wellList, std::vector<std::string>& in_categoryList,std::unordered_map<std::string,std::string>& in_refseqToGene, std::unordered_map<std::string,uint32_t>& in_erccToIndex,
		std::unordered_map<std::string,uint32_t>& in_geneToIndex ) :

		markMultiHits(in_markMultiHits),
		markUnknowns(in_markUnknowns),
		mixtureOfWells(in_mixtureOfWells),
		properReads(in_properReads),
		barcodeSize(in_barcodeSize),
		umiSize(in_umiSize),
		nthreads(barcodeSize ? in_nthreads : 1),
		erccListSize(maxErccListSize),
		geneListSize(maxGeneListSize),
		nwells (barcodeSize ? NWELLS : 1),
		ncategories(maxGeneListSize + maxErccListSize + 1),
		alignsDir(in_alignsDir),
		sampleId(in_sampleId),
		wellList(in_wellList),
	    nbinbits(in_nbinbits),
	    binsizebits(in_binsizebits),
		leftBitMask(1 << sizeof(T)-1),
		umibits(2*umiSize),
		posMask((1 << nbinbits) -1),
		categoryList(in_categoryList),
		refseqToGene(in_refseqToGene),
		erccToIndex(in_erccToIndex),
		geneToIndex(in_geneToIndex)	
		{
			//non-gene dependent counts
		    total_reads = new uint32_t [nthreads];
		    assigned_reads = new uint32_t [nthreads];
			aligned_reads = maxGeneListSize ? new uint32_t [nthreads]: 0;
			assigned_unknown_reads = markUnknowns ? new uint32_t [nthreads] : 0;
			assigned_unknown_umi = (markUnknowns && umiSize) ? new uint32_t [nthreads] : 0;
			if (markMultiHits){
			    total_reads_mm = new uint32_t [nthreads];
			    assigned_reads_mm = new uint32_t [nthreads];
				aligned_reads_mm = maxGeneListSize ? new uint32_t [nthreads] : 0;
				assigned_unknown_reads_mm = markUnknowns ? new uint32_t [nthreads] : 0;
				assigned_unknown_umi_mm = (markUnknowns && umiSize) ? new uint32_t [nthreads] : 0;
			}
			else{
			    total_reads_mm = 0;
			    assigned_reads_mm = 0;
				aligned_reads_mm = 0;
				assigned_unknown_reads_mm = 0;
				assigned_unknown_umi_mm =  0;				
			}
			
			uint32_t *nonCategorySums[] = {total_reads, assigned_reads, aligned_reads, assigned_unknown_reads, assigned_unknown_umi,total_reads_mm, assigned_reads_mm, aligned_reads_mm, assigned_unknown_reads_mm, assigned_unknown_umi_mm};
	
			for (uint32_t *nonCategorySum : nonCategorySums){
				if (nonCategorySum) memset(nonCategorySum,0,nthreads*sizeof(uint32_t));
			} 

			//allocate memory for gene/category based counts
			//first dimension is nwells - stores pointers to block of nwells*ncategories
			categoryTotal = new uint32_t* [nwells];
			categoryUmi = umiSize? new uint32_t* [nwells] : 0;
			unknown_set_w = markUnknowns ? new std::unordered_set<std::string> [nwells]() : 0;
			unknown_set_w_mm = (markUnknowns && markMultiHits) ? new std::unordered_set<std::string> [nwells]() : 0;
			
			if (markMultiHits){
				categoryTotal_mm = new uint32_t* [nwells];
				categoryUmi_mm = umiSize? new uint32_t* [nwells] : 0;
			}
			else{
				categoryTotal_mm = 0;
				categoryUmi_mm =  0;					
			}
			uint32_t **categorySums[]={categoryTotal,categoryUmi,categoryTotal_mm,categoryUmi_mm};
			for (uint32_t **categorySum : categorySums){
				if (categorySum){
					categorySum[0] = new uint32_t[ncategories*nwells];
					memset(categorySum[0],0,ncategories*nwells*sizeof(uint32_t));
					for(int i=1;i<nwells;i++){	   
						categorySum[i]=categorySum[i-1]+ncategories;
					}
				}
			}

		}
	
		void print(std::string resultsDir,bool filterSAMoutput){
			if (filterSAMoutput) merge_saf();
			else merge_sam();
			
			std::string logFile=resultsDir+"/"+sampleId+".unq.log.dat";
			std::string unknownFile=resultsDir+"/"+sampleId+".unq.unknown_list";
			std::string geneFile=resultsDir+"/"+sampleId+".unq.refseq.total.dat";
			std::string geneUmiFile=resultsDir+"/"+sampleId+".unq.refseq.umi.dat";
			std::string spikeFile=resultsDir+"/"+sampleId+".unq.spike.total.dat";
			std::string spikeUmiFile=resultsDir+"/"+sampleId+".unq.spike.umi.dat";
			std::string wellTotalFile=resultsDir+"/"+sampleId+".unq.well_summary.dat";
			std::string logFile_mm=resultsDir+"/"+sampleId+".all.log.dat";
			std::string unknownFile_mm=resultsDir+"/"+sampleId+".all.unknown_list";
			std::string geneFile_mm=resultsDir+"/"+sampleId+".all.refseq.total.dat";
			std::string geneUmiFile_mm=resultsDir+"/"+sampleId+".all.refseq.umi.dat";
			std::string spikeFile_mm=resultsDir+"/"+sampleId+".all.spike.total.dat";
			std::string spikeUmiFile_mm=resultsDir+"/"+sampleId+".all.spike.umi.dat";
			std::string wellTotalFile_mm=resultsDir+"/"+sampleId+".all.well_summary.dat";
		 
		  //combine unknown sets - assume that bad wells are ignored
			if (markUnknowns){
				for(int w=0; w< nwells; w++){
					for (auto itr = unknown_set_w[w].begin(); itr != unknown_set_w[w].end(); ++itr){
					 unknown_set.insert(*itr);
					}
				}
				for(int w=0; w< nwells; w++){
					for (auto itr = unknown_set_w_mm[w].begin(); itr != unknown_set_w_mm[w].end(); ++itr){
						unknown_set_mm.insert(*itr);
					}
				}
			}
	
			FILE *fp=fopen(logFile.c_str(),"w");
			if(!fp)exit(EXIT_FAILURE);
			fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
			fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				sampleId.c_str(),
				sumThreadCounts(total_reads),
				sumThreadCounts(assigned_reads),
				sumThreadCounts(aligned_reads),
				sumCategoryCounts(categoryTotal,geneListSize,erccListSize),
				sumCategoryCounts(categoryUmi,geneListSize,erccListSize),
				sumCategoryCounts(categoryTotal,geneListSize+erccListSize,1),
				sumCategoryCounts(categoryUmi,geneListSize+erccListSize,1),										
				sumCategoryCounts(categoryTotal,0,geneListSize),
				sumCategoryCounts(categoryUmi,0,geneListSize),
				sumThreadCounts(assigned_unknown_reads),
				sumThreadCounts(assigned_unknown_umi)
			);
			fclose(fp);
			if (markMultiHits){
				if(!(fp=fopen(logFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
				fprintf(fp,		"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
				fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					sampleId.c_str(),
					sumThreadCounts(total_reads_mm),
					sumThreadCounts(assigned_reads_mm),
					sumThreadCounts(aligned_reads_mm),
					sumCategoryCounts(categoryTotal_mm,geneListSize,erccListSize),
					sumCategoryCounts(categoryUmi_mm,geneListSize,erccListSize),
					sumCategoryCounts(categoryTotal_mm,geneListSize+erccListSize,1),
					sumCategoryCounts(categoryUmi_mm,geneListSize+erccListSize,1),										
					sumCategoryCounts(categoryTotal_mm,0,geneListSize),
					sumCategoryCounts(categoryUmi_mm,0,geneListSize),
					sumThreadCounts(assigned_unknown_reads_mm),
					sumThreadCounts(assigned_unknown_umi_mm)
				);
				fclose(fp);
			}
	
			if (markUnknowns){
				if(!(fp=fopen(unknownFile.c_str(),"w")))exit(EXIT_FAILURE);
				for (auto itr = unknown_set.begin(); itr != unknown_set.end(); ++itr) {
					fprintf(fp,"%s\n",itr->c_str());
				}
				fclose(fp);
				if (markMultiHits){
					if(!(fp=fopen(unknownFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
					for (auto itr = unknown_set_mm.begin(); itr != unknown_set_mm.end(); ++itr) {
						fprintf(fp,"%s\n",itr->c_str());
					}
					fclose(fp);
				}
			}
			printCounts(geneFile,0,geneListSize,categoryTotal);
			printCounts(geneFile_mm,0,geneListSize,categoryTotal_mm);
			printCounts(geneUmiFile,0,geneListSize,categoryUmi);
			printCounts(geneUmiFile_mm,0,geneListSize,categoryUmi_mm);
		 
			printCounts(spikeFile,geneListSize,erccListSize,categoryTotal);
			printCounts(spikeFile_mm,geneListSize,erccListSize,categoryTotal_mm);
			printCounts(spikeUmiFile,geneListSize,erccListSize,categoryUmi);
			printCounts(spikeUmiFile_mm,geneListSize,erccListSize,categoryUmi_mm);
	 
			//print well totals
			if(!(fp=fopen(wellTotalFile.c_str(),"w")))exit(EXIT_FAILURE);
			for(int i=0;i<nwells;i++){
				fprintf(fp,"\t%s",wellList[i].c_str());
			}
			fprintf(fp,"\n");
			printWellSums(fp,"Refseq_Total",categoryTotal,0,geneListSize);
			printWellSums(fp,"Refseq_Umi",categoryUmi,0,geneListSize);
			printWellSums(fp,"Spike_Total",categoryTotal,geneListSize,erccListSize);
			printWellSums(fp,"Spike_Umi",categoryUmi,geneListSize,erccListSize);
			fclose(fp);
			if(markMultiHits){
				if(!(fp=fopen(wellTotalFile_mm.c_str(),"w")))exit(EXIT_FAILURE); 
				for(int i=0;i<wellList.size();i++){
					fprintf(fp,"\t%s",wellList[i].c_str());
				}
				fprintf(fp,"\n");
				printWellSums(fp,"Refseq_Total",categoryTotal_mm,0,geneListSize);
				printWellSums(fp,"Refseq_Umi",categoryUmi_mm,0,geneListSize);
				printWellSums(fp,"Spike_Total",categoryTotal_mm,geneListSize,erccListSize);
				printWellSums(fp,"Spike_Umi",categoryUmi_mm,geneListSize,erccListSize);
				fclose(fp);
			}		
		}
	

		~Counts(){
			//clean up 
			delete[] categoryTotal[0];
			delete[] categoryUmi[0];   
			delete[] categoryTotal;
			delete[] categoryUmi;   
			delete[] total_reads;
			delete[] assigned_reads;
			delete[] aligned_reads;	
			delete[] assigned_unknown_reads;
			delete[] assigned_unknown_umi;
			delete[] total_reads_mm;
			delete[] assigned_reads_mm;
			delete[] aligned_reads_mm;	
			delete[] assigned_unknown_reads_mm;
			delete[] assigned_unknown_umi_mm;
		}

};	

 
int main(int argc, char *argv[]){
	
	bool markMultiHits=0,markUnknowns=0,properReads=0;
	uint8_t barcodeSize=6, umiSize=10, nbinbits=16, binsizebits=0;
	bool geneLevelFilter=0,filteredSAMfiles=0,mixtureOfWells=0;
	std::string sampleId="",sym2ref="", ercc_fasta="", barcodes="", alignDir="", resultsDir="",countsFile="";
	int opt,verbose=0,nthreads=1;
	struct optparse options;
	optparse_init(&options, argv);	
 
	std::string errmsg="umimerge_parallel v?hfmn:go:t:i:s:e:m:c:b:a:p:\n-h -?  (display this message)\n-v (Verbose mode)\n-f Merge filtered SAM files (*.saf) instead of full SAM files\n-s <sym2ref file>\n-e <ercc_fasta file>\n-b <barcode_file>\n-a <aligns directory>\n-o <resultsDir(counts directory)>\n-c <binary file with umicounts - default is UMIcounts.bin in the input SAM directory>\n-t <number of threads (1)>\n-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... \n-m The aligned files in the output directory are not separated into wells by subdirectories and the filename is used to identify the well\n\nRequired params are -i sampleId -s sym2ref -e ercc_fasta -b barcodes -a alignDir -o resultsDir\n\nExample:\n\numimerge_parallel -i RNAseq_20150409 -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -a Aligns -o Counts -t 4 -p 0\n";
	while ((opt = optparse(&options, "v?hfmn:go:t:i:s:e:m:c:b:a:p:")) != -1) {
		switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'n':
			 nbinbits=atoi(options.optarg);
			 posMask=(1 << nbinbits) -1;
			 //set nbinbits - default is 16 bits otherwise  
			break;
			case 'f':
			 //use filtered SAM files and not sam files
			 filteredSAMfiles=1;
			break;
			case 'g':
				//set gene level filtering
				nbinbits=0;
			break;
			case 'o':
			 resultsDir=std::string(options.optarg);
			break;
			case 't':
				nthreads=atoi(options.optarg);
			break;
			case 'i':
				sampleId=std::string(options.optarg);
			break;  
			case 's':
				sym2ref=std::string(options.optarg);
			break;
			case 'e':
				ercc_fasta=std::string(options.optarg);
			break; 
			case 'b':
				barcodes=std::string(options.optarg);
			break;
			case 'm':
				mixtureOfWells=1;
			break;   
			case 'a':
				alignDir=std::string(options.optarg);
			break;
			case 'p':		 
				binsizebits=atoi(options.optarg);
			break;
			case 'P':		 
				properReads=1;
				//for paired ends 
			break;
			case  'M':
				markMultiHits=1;
			break;
			case  'u':
				markUnknowns=1;
			break;
			case  'B':
				barcodeSize=atoi(options.optarg);;
			break;
			case  'U':
				umiSize=atoi(options.optarg);
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
	if(sampleId =="" || sym2ref=="" ||  alignDir=="" || resultsDir==""){
		fprintf(stderr,"Required params are -i sampleId -s sym2ref  -a alignDir -o resultsDir\n");
		exit(EXIT_FAILURE);
	}	
 
	std::unordered_map<std::string,std::string> refseqToGene;
	std::vector<std::string>erccList, geneList,wellList,categoryList;
	std::unordered_map<std::string,uint32_t>wellToIndex;
	std::unordered_map<std::string,uint32_t>erccToIndex;
	std::unordered_map<std::string,uint32_t>geneToIndex;
	std::vector<std::string> unknown_list;
	
	if (barcodeSize) readWells(barcodes,wellToIndex,wellList);
	else wellList.push_back("Total");
	
	readERCC(ercc_fasta,erccList);
	readRefseq(sym2ref,refseqToGene ,geneList);
 
	categoryList=geneList;
	if (erccList.size()) categoryList.insert(categoryList.end(),erccList.begin(), erccList.end());
	categoryList.push_back("chrM");
 
	for(int i=0;i<geneList.size();i++)
		geneToIndex[geneList[i]]=i;
	for(int i=0;i<erccList.size();i++)
		erccToIndex[erccList[i]]=i;
	
	uint8_t bitSize=bitsNeeded(categoryList.size(),umiSize, nbinbits, markMultiHits);
	//outputfile names
	
    if (bitSize > 64){
		fprintf(stderr,"Maximum hash size of 64 bits exceeded - %d bits required for encoding\n", bitSize);
		exit(EXIT_FAILURE);		
	}
	if (bitSize > 32){
		Counts <uint64_t> count (nthreads, erccList.size(),geneList.size(), markMultiHits,mixtureOfWells, markUnknowns,properReads,barcodeSize, umiSize, nbinbits, binsizebits, alignDir, sampleId, wellList, categoryList, refseqToGene, erccToIndex, geneToIndex );
		count.print(resultsDir,filterSAMoutput);	
    }
    else if (bitSize > 16){
		Counts <uint32_t> count (nthreads, erccList.size(),geneList.size(), markMultiHits,mixtureOfWells, markUnknowns,properReads,barcodeSize, umiSize, nbinbits, binsizebits, alignDir, sampleId, wellList, categoryList, refseqToGene, erccToIndex, geneToIndex );
		count.print(resultsDir,filterSAMoutput);	
	}
    else {
		Counts <uint16_t> count (nthreads, erccList.size(),geneList.size(), markMultiHits,mixtureOfWells, markUnknowns,properReads,barcodeSize, umiSize, nbinbits, binsizebits, alignDir, sampleId, wellList, categoryList, refseqToGene, erccToIndex, geneToIndex );
		count.print(resultsDir,filterSAMoutput);		
	}	 
}
