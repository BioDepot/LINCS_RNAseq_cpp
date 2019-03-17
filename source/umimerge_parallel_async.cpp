#include <zlib.h>  
#include <omp.h>
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

//asynch option will take a lot more memory as all hashes have to be formed and kept in memory for UMI hashing
//OK for 96 wells but would crash and burn for 384 wells

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
bool doneFileExists(const std::string &doneFile);

class Counts{
	private:
		vector<unordered_set<uint64_t>> umi_seens(NWELLS+1);
		vector<unordered_set<uint64_t>> umi_seens_mm(NWELLS+1);
	public:
		unsigned int
		*total_reads,
		*assigned_reads,
		*assigned_aligned_reads,
		*assigned_mito_reads,
		*assigned_mito_umi,
		*assigned_unknown_reads,
		*assigned_unknown_umi,
		*total_reads_mm,
		*assigned_reads_mm,
		*assigned_aligned_reads_mm,	
		*assigned_mito_reads_mm,
		*assigned_mito_umi_mm,
		*assigned_unknown_reads_mm,
		*assigned_unknown_umi_mm ,
		*spike_total[NWELLS],
		*spike_umi[NWELLS],   
		*total[NWELLS],
		*umi[NWELLS],   
		*spike_total_mm[NWELLS],
		*spike_umi_mm[NWELLS],   
		*total_mm[NWELLS],
		*umi_mm[NWELLS];

		
		//this hash is necessary to join the unknown hashes	 
		unordered_set<string>unknown_set;
		unordered_set<string>unknown_set_w[NWELLS];
		unordered_set<string>unknown_set_mm;
		unordered_set<string>unknown_set_mm_w[NWELLS];

		Counts(unsigned int nThreads,unsigned int erccListSize, unsigned int geneListSize){
			init_counts(nThreads,erccListSize,geneListSize);
		}	
	 
		void init_counts(unsigned int nThreads,unsigned int erccListSize, unsigned int geneListSize){

			//allocate memory for total
			spike_total[0]=new unsigned int[erccListSize*NWELLS];
			spike_umi[0]=new unsigned int[erccListSize*NWELLS];      
			spike_total_mm[0]=new unsigned int[erccListSize*NWELLS];
			spike_umi_mm[0]=new unsigned int[erccListSize*NWELLS];   
		
			total[0]=new unsigned int[geneListSize*NWELLS];
			umi[0]=new unsigned int[geneListSize*NWELLS];
			total_mm[0]=new unsigned int[geneListSize*NWELLS];
			umi_mm[0]=new unsigned int[geneListSize*NWELLS];
	

			//unq counts
			total_reads = new unsigned int [nThreads],
			assigned_reads = new unsigned int [nThreads],
			assigned_aligned_reads = new unsigned int [nThreads],	
			assigned_mito_reads = new unsigned int [nThreads],
			assigned_mito_umi = new unsigned int [nThreads],
			assigned_unknown_reads = new unsigned int [nThreads],
			assigned_unknown_umi = new unsigned int [nThreads],
  
			//mm *counts
			total_reads_mm = new unsigned int [nThreads],
			assigned_reads_mm = new unsigned int [nThreads],
			assigned_aligned_reads_mm = new unsigned int [nThreads],	
			assigned_mito_reads_mm = new unsigned int [nThreads],
			assigned_mito_umi_mm = new unsigned int [nThreads],
			assigned_unknown_reads_mm = new unsigned int [nThreads],
			assigned_unknown_umi_mm = new unsigned int [nThreads],
  
			//2D arrays for spike ins
			//arrays for pointers
			//allocate memory for total
			spike_total[0]=new unsigned int[erccListSize*NWELLS];
			spike_umi[0]=new unsigned int[erccListSize*NWELLS];      
			spike_total_mm[0]=new unsigned int[erccListSize*NWELLS];
			spike_umi_mm[0]=new unsigned int[erccListSize*NWELLS];   
		
			total[0]=new unsigned int[geneListSize*NWELLS];
			umi[0]=new unsigned int[geneListSize*NWELLS];
			total_mm[0]=new unsigned int[geneListSize*NWELLS];
			umi_mm[0]=new unsigned int[geneListSize*NWELLS];
		
			//set to zero
		
			memset(spike_total[0],0,erccListSize*NWELLS*sizeof(unsigned int));
			memset(spike_umi[0],0,erccListSize*NWELLS*sizeof(unsigned int));      
			memset(spike_total_mm[0],0,erccListSize*NWELLS*sizeof(unsigned int));
			memset(spike_umi_mm[0],0,erccListSize*NWELLS*sizeof(unsigned int));   
		
			memset(total[0],0,geneListSize*NWELLS*sizeof(unsigned int));
			memset(umi[0],0,geneListSize*NWELLS*sizeof(unsigned int));
			memset(total_mm[0],0,geneListSize*NWELLS*sizeof(unsigned int));
			memset(umi_mm[0],0,geneListSize*NWELLS*sizeof(unsigned int));
		
			memset(total_reads,0,nThreads*sizeof(unsigned int));
			memset(assigned_reads,0,nThreads*sizeof(unsigned int));
			memset(assigned_aligned_reads,0,nThreads*sizeof(unsigned int));	
			memset(assigned_mito_reads,0,nThreads*sizeof(unsigned int));
			memset(assigned_mito_umi,0,nThreads*sizeof(unsigned int));
			memset(assigned_unknown_reads,0,nThreads*sizeof(unsigned int));
			memset(assigned_unknown_umi,0,nThreads*sizeof(unsigned int));
			memset(total_reads_mm,0,nThreads*sizeof(unsigned int));
			memset(assigned_reads_mm,0,nThreads*sizeof(unsigned int));
			memset(assigned_aligned_reads_mm,0,nThreads*sizeof(unsigned int));	
			memset(assigned_mito_reads_mm,0,nThreads*sizeof(unsigned int));
			memset(assigned_mito_umi_mm,0,nThreads*sizeof(unsigned int));
			memset(assigned_unknown_reads_mm,0,nThreads*sizeof(unsigned int));
			memset(assigned_unknown_umi_mm,0,nThreads*sizeof(unsigned int));	
  
			//initialize pointers
			for(int i=1;i<NWELLS;i++){
				spike_total[i]=spike_total[i-1]+erccListSize;
				spike_umi[i]=spike_umi[i-1]+erccListSize;      
				spike_total_mm[i]=spike_total_mm[i-1]+erccListSize;
				spike_umi_mm[i]=spike_umi_mm[i-1]+erccListSize;   
				total[i]=total[i-1]+geneListSize;
				umi[i]=umi[i-1]+geneListSize;
				total_mm[i]=total_mm[i-1]+geneListSize;
				umi_mm[i]=umi_mm[i-1]+geneListSize;
			}
		}
		
		void print(string dge_dir, string sample_id, unsigned int nThreads, vector<string> &erccList ,vector<string> &geneList, vector<string> &wellList){
			string umiBarcodesFile=dge_dir+"/"+sample_id+".unq.umi.dat";
			string logFile=dge_dir+"/"+sample_id+".unq.log.dat";
			string unknownFile=dge_dir+"/"+sample_id+".unq.unknown_list";
			string totalAlignFile=dge_dir+"/"+sample_id+".unq.refseq.total.dat";
			string umiFile=dge_dir+"/"+sample_id+".unq.refseq.umi.dat";
			string erccFile=dge_dir+"/"+sample_id+".unq.spike.total.dat";
			string erccUMIFile=dge_dir+"/"+sample_id+".unq.spike.umi.dat";
			string wellTotalFile=dge_dir+"/"+sample_id+".unq.well_summary.dat";
	 
			string logFile_mm=dge_dir+"/"+sample_id+".all.log.dat";
			string unknownFile_mm=dge_dir+"/"+sample_id+".all.unknown_list";
			string totalAlignFile_mm=dge_dir+"/"+sample_id+".all.refseq.total.dat";
			string umiFile_mm=dge_dir+"/"+sample_id+".all.refseq.umi.dat";
			string erccFile_mm=dge_dir+"/"+sample_id+".all.spike.total.dat";
			string erccUMIFile_mm=dge_dir+"/"+sample_id+".all.spike.umi.dat";
			string wellTotalFile_mm=dge_dir+"/"+sample_id+".all.well_summary.dat";
	 
			//combine unknown sets - assume that bad wells are ignored
			for(int w=0; w< NWELLS; w++){
				for (auto itr = unknown_set_w[w].begin(); itr != unknown_set_w[w].end(); ++itr){
					unknown_set.insert(*itr);
				}
			}
			for(int w=0; w< NWELLS; w++){
				for (auto itr = unknown_set_mm_w[w].begin(); itr != unknown_set_mm_w[w].end(); ++itr){
					unknown_set_mm.insert(*itr);
				}
			}
			FILE *fp=fopen(logFile.c_str(),"w");
			if(!fp)exit(EXIT_FAILURE);
			fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
			fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					sample_id.c_str(),
					sumCounts(total_reads,nThreads),
					sumCounts(assigned_reads,nThreads),
					sumCounts(assigned_aligned_reads,nThreads),
					sumCounts(spike_total,NWELLS,erccList.size()),
					sumCounts(spike_umi,NWELLS,erccList.size()),
					sumCounts(assigned_mito_reads,nThreads),
					sumCounts(assigned_mito_umi,nThreads),
					sumCounts(total,NWELLS,geneList.size()),
					sumCounts(umi,NWELLS,geneList.size()),
					sumCounts(assigned_unknown_reads,nThreads),
					sumCounts(assigned_unknown_umi,nThreads)
			);
			fclose(fp);
			if(!(fp=fopen(logFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
			fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
			fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				sample_id.c_str(),
				sumCounts(total_reads_mm,nThreads),
				sumCounts(assigned_reads_mm,nThreads),
				sumCounts(assigned_aligned_reads_mm,nThreads),		
				sumCounts(spike_total_mm,NWELLS,erccList.size()),
				sumCounts(spike_umi_mm,NWELLS,erccList.size()),
				sumCounts(assigned_mito_reads_mm,nThreads),
				sumCounts(assigned_mito_umi_mm,nThreads),
				sumCounts(total_mm,NWELLS,geneList.size()),
				sumCounts(umi_mm,NWELLS,geneList.size()),
				sumCounts(assigned_unknown_reads_mm,nThreads),
				sumCounts(assigned_unknown_umi_mm,nThreads)
			);
			fclose(fp);
	 
			if(!(fp=fopen(unknownFile.c_str(),"w")))exit(EXIT_FAILURE);
			for (auto itr = unknown_set.begin(); itr != unknown_set.end(); ++itr) {
				fprintf(fp,"%s\n",itr->c_str());
			}
			fclose(fp);
			if(!(fp=fopen(unknownFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
			for (auto itr = unknown_set_mm.begin(); itr != unknown_set_mm.end(); ++itr) {
				fprintf(fp,"%s\n",itr->c_str());
			}
			fclose(fp);
	 
			if(!(fp=fopen(totalAlignFile.c_str(),"w")))exit(EXIT_FAILURE);
			for(int i=0;i<wellList.size();i++){
				fprintf(fp,"\t%s",wellList[i].c_str());
			}
			fprintf(fp,"\n");
			for(int i=0;i<geneList.size();i++){
				fprintf(fp,"%s\t",geneList[i].c_str());
				for(int j=0;j<NWELLS-1;j++){
					fprintf(fp,"%d\t",total[j][i]);
				}
				fprintf(fp,"%d\n",total[NWELLS-1][i]);
			} 
			fclose(fp);	
			if(!(fp=fopen(totalAlignFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
			for(int i=0;i<wellList.size();i++){
				fprintf(fp,"\t%s",wellList[i].c_str());
			}
			fprintf(fp,"\n");
			for(int i=0;i<geneList.size();i++){
				fprintf(fp,"%s\t",geneList[i].c_str());
				for(int j=0;j<NWELLS-1;j++){
					fprintf(fp,"%d\t",total_mm[j][i]);
				}
				fprintf(fp,"%d\n",total_mm[NWELLS-1][i]);
			} 
			fclose(fp);
	 
			if(!(fp=fopen(umiFile.c_str(),"w")))exit(EXIT_FAILURE);
			for(int i=0;i<wellList.size();i++){
				fprintf(fp,"\t%s",wellList[i].c_str());
			}
			fprintf(fp,"\n");  
			for(int i=0;i<geneList.size();i++){
				fprintf(fp,"%s\t",geneList[i].c_str());
			for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",umi[j][i]);
			}
			fprintf(fp,"%d\n",umi[NWELLS-1][i]);
		}  
		fclose(fp);
		if(!(fp=fopen(umiFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n"); 
		for(int i=0;i<geneList.size();i++){
			fprintf(fp,"%s\t",geneList[i].c_str());
			for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",umi_mm[j][i]);
			}
			fprintf(fp,"%d\n",umi_mm[NWELLS-1][i]);
		}
		fclose(fp);
		if(!(fp=fopen(erccFile.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n"); 
		for(int i=0;i<erccList.size();i++){
			fprintf(fp,"%s",erccList[i].c_str());
			for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",spike_umi[j][i]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		if(!(fp=fopen(erccFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n"); 
		for(int i=0;i<erccList.size();i++){
			fprintf(fp,"%s",erccList[i].c_str());
			for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",spike_umi_mm[j][i]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);   
	 
		if(!(fp=fopen(erccUMIFile.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
		for(int i=0;i<erccList.size();i++){
			fprintf(fp,"%s",erccList[i].c_str());
			for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",spike_total[j][i]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		if(!(fp=fopen(erccUMIFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n"); 
		for(int i=0;i<erccList.size();i++){
			fprintf(fp,"%s",erccList[i].c_str());
			for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",spike_total_mm[j][i]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		if(!(fp=fopen(wellTotalFile.c_str(),"w")))exit(EXIT_FAILURE);
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 
		fprintf(fp,"Refseq_Total"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(total,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Refseq_UMI"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(umi,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Spike_Total"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(spike_total,i,erccList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Spike_UMI"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(spike_umi,i,erccList.size()));
		}
		fprintf(fp,"\n");
		fclose(fp);
	 
		if(!(fp=fopen(wellTotalFile_mm.c_str(),"w")))exit(EXIT_FAILURE); 
		for(int i=0;i<wellList.size();i++){
			fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 
		fprintf(fp,"Refseq_Total"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(total_mm,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Refseq_UMI"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(umi_mm,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Spike_Total"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(spike_total_mm,i,erccList.size()));
		}
		fprintf(fp,"\n");
		
		fprintf(fp,"Spike_UMI"); 
		for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(spike_umi_mm,i,erccList.size()));
		}
		fprintf(fp,"\n");
		fclose(fp);		
	}
 
	template <class T1,class T2> void merge_parallel (int nThreads,vector<string> &erccList ,vector<string> &geneList, vector<string> &wellList, string aligned_dir, umipanel<T1,T2> &barcodePanel,  unordered_map<string,string> refseq_to_gene,unordered_map<string,unsigned int>well_to_index,unordered_map<string,unsigned int>ercc_to_index,unordered_map<string,unsigned int>gene_to_index, uint32_t posMask, int binSize, int nbins, bool geneLevelFilter){
													  const uint64_t nCategories=geneList.size()+erccList.size()+1; //add1 for chrM
	   //fprintf (stderr, "nCategories is %d Gene list size is %d errccList.szie is %d \n",nCategories, geneList.size(),erccList.size());
		#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
		for	(int wellIndex=0; wellIndex<NWELLS; wellIndex++){
				
		//hashes
		unordered_set<uint64_t>umi_seen;
		unordered_set<uint64_t>umi_seen_mm;
		unordered_set<string>unknown_umi_seen;
		unordered_set<string>unknown_umi_seen_mm;
	
	 	vector<string> inputFiles;
	 	glob_t glob_result;
	
		string well=barcodePanel.wells[wellIndex];
		string globString=aligned_dir+"/"+well+"/"+"*.sam";
		glob(globString.c_str(),GLOB_TILDE,NULL,&glob_result);
		for(int j=0; j<glob_result.gl_pathc; j++){
			inputFiles.push_back(string(glob_result.gl_pathv[j]));
		}
		for(int i=0;i<inputFiles.size();i++){
			const unsigned int tid=omp_get_thread_num();
			string fullLine;
			fprintf(stderr,"thread %d of %d working on %s\n",tid,nThreads,inputFiles[i].c_str());
			std::ifstream fstream(inputFiles[i]);
			if(!fstream.is_open()){		fprintf(stderr,"unable to open file %s\n", inputFiles[i].c_str());continue;}
			while(getline(fstream,fullLine)){
				 if(fullLine[0] != '@'){
					total_reads[tid]++;
					total_reads_mm[tid]++;;
		   
				   //get barcode
					vector <string> items;
					vector <string> tempItems;
					splitStr(fullLine," \t",items);
	
					splitStr(items[0],":",tempItems);
					string fullBarcode=tempItems[tempItems.size()-1];
			 	  
					//string well=tempItems[tempItems.size()-2];
			
					string barcode=fullBarcode.substr(6,UMISIZE); //change if barcode size changes
				   
					//ambigCheck
					//skip if ambiguous barcode and get unique barcode index from sequence
					unsigned int barcodeIndex=0;
					if(ambigCheck(barcode,barcodeIndex))continue;
					assigned_reads[tid]++;
					assigned_reads_mm[tid]++;
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
						assigned_aligned_reads[tid]++;
						assigned_aligned_reads_mm[tid]++;
							
						if (aligned_id.substr(0,4) == "ERCC"){
							const int erccIndex=ercc_to_index[aligned_id];
							uint64_t umicode=barcodeIndex*nCategories+erccIndex+geneList.size();
							if(!geneLevelFilter){
								uint32_t positionCode=(pos >> binSize) & posMask; 
								umicode = (umicode << nbins) | positionCode;
							}
							if(!best_hits_list.size()){ 
								spike_total[wellIndex][erccIndex]++;						
								if(!umi_seen.count(umicode)){
									umi_seen.insert(umicode);
									spike_umi[wellIndex][erccIndex]++;
								}
							} 
							spike_total_mm[wellIndex][erccIndex]++;
							if(!umi_seen_mm.count(umicode)){
								umi_seen_mm.insert(umicode);
								spike_umi_mm[wellIndex][erccIndex]++;
							}
						}
						else if 	(aligned_id.substr(0,4) == "chrM"){
							uint64_t umicode=barcodeIndex*nCategories+nCategories-1;
							if(!geneLevelFilter){
								uint32_t positionCode=(pos >> binSize) & posMask; 
								umicode = (umicode << nbins) | positionCode;
							}
							if(!best_hits_list.size()){
								assigned_mito_reads[tid]++;
								if(!umi_seen.count(umicode)){
								 	assigned_mito_umi[tid]++;
								 	umi_seen.insert(umicode);
								}	
							} 
							assigned_mito_reads_mm[tid]++;
							if(!umi_seen_mm.count(umicode)){
								assigned_mito_umi_mm[tid]++;
								umi_seen_mm.insert(umicode);
							}
						}
						else if (refseq_to_gene.count(aligned_id)){
							string gene = refseq_to_gene[aligned_id];
							const unsigned int geneIndex =gene_to_index[gene];
							uint64_t umicode=barcodeIndex*nCategories+geneIndex;
							if(!geneLevelFilter){
								uint32_t positionCode=(pos >> binSize) & posMask; 
								umicode = (umicode << nbins) | positionCode;
							}
							bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
							if(!multiGene){
								total[wellIndex][geneIndex]++;
								if(!umi_seen.count(umicode)){
									umi[wellIndex][geneIndex]++;
									umi_seen.insert(umicode);
								}
							}
							total_mm[wellIndex][geneIndex]++;
							if(!umi_seen_mm.count(umicode)){
								umi_mm[wellIndex][geneIndex]++;
								umi_seen_mm.insert(umicode);
							}					
						}
						else{
							string umicode=barcodePanel.wells[wellIndex] + barcode + aligned_id;
							if(!best_hits_list.size()){
								assigned_unknown_reads[tid]++;
								if(!unknown_set_w[wellIndex].count(aligned_id)){
								 	unknown_set.insert(aligned_id);
								}
								if(!unknown_umi_seen.count(umicode)){
									assigned_unknown_umi[tid]++;	
									unknown_umi_seen.insert(umicode);
								}
							}
							if(!unknown_set_mm.count(aligned_id)){
								unknown_set_mm_w[wellIndex].insert(aligned_id);
							}
							assigned_unknown_reads_mm[tid]++;
							if(!unknown_umi_seen_mm.count(umicode)){
								assigned_unknown_umi_mm[tid]++;
								unknown_umi_seen_mm.insert(umicode);
							}				
						}		 
					}
				}
				fstream.close();
			}
		}
	}
	
 
	template <class T1,class T2> void merge_filter(int nThreads,vector<string> &erccList ,vector<string> &geneList, vector<string> &wellList, string aligned_dir, umipanel<T1,T2> &barcodePanel,  unordered_map<string,string> refseq_to_gene,unordered_map<string,unsigned int>well_to_index,unordered_map<string,unsigned int>ercc_to_index,unordered_map<string,unsigned int>gene_to_index, uint32_t posMask, int binSize, int nbins, bool geneLevelFilter){
								const uint64_t nCategories=geneList.size()+erccList.size()+1; //add1 for chrM
		
		//each thread works on a well


		#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
		for	(int wellIndex=0; wellIndex<NWELLS; wellIndex++){
			uint64_t mask63 = (1ULL << 63) -1; //all 1's except for left most bit
		  
		 	//find matching files in the separated wells directory - we skip X - i.e. unmappables now
	 		vector<string> inputFiles;
			glob_t glob_result;
			 //string well=(wellIndex >=NWELLS) ? "X" :barcodePanel.wells[wellIndex];
			string well=barcodePanel.wells[wellIndex];
			string globString=aligned_dir+"/"+well+"/"+"*.saf";
			glob(globString.c_str(),GLOB_TILDE,NULL,&glob_result);
			for(int j=0; j<glob_result.gl_pathc; j++){
				inputFiles.push_back(string(glob_result.gl_pathv[j]));
			}
			for(int i=0;i<inputFiles.size();i++){
				const unsigned int tid=omp_get_thread_num();
				uint64_t umicode;
				fprintf(stderr,"thread %d of %d working on %s\n",tid,nThreads,inputFiles[i].c_str());
				//read in each file
				FILE *fp=fopen(inputFiles[i].c_str(),"r");
				while (fread(&umicode,sizeof(uint64_t),1,fp)){
					uint32_t category;
					bool multiGeneHit=(bool) (umicode >> 63);
					umicode=umicode & mask63;
					if(geneLevelFilter) category=umicode % nCategories;
					else category = (umicode >> nbins) % nCategories;
					assigned_aligned_reads[tid]++;
					assigned_aligned_reads_mm[tid]++;
						
					if (category == nCategories -1){
							//mitochondrial
						assigned_mito_reads_mm[tid]++;
						if (!multiGeneHit){
							assigned_mito_reads[tid]++;
							if(!umi_seen.count(umicode)){
								umi_seen.insert(umicode);
								assigned_mito_umi[tid]++;
							}
						}
						if(!umi_seen_mm.count(umicode)){
							umi_seen_mm.insert(umicode);
							assigned_mito_umi_mm[tid]++;
						}
					}
					else if (category >= geneList.size()){
						//spike
						const unsigned int erccIndex=category-geneList.size();
						spike_total_mm[wellIndex][erccIndex]++;
						if (!multiGeneHit){
							spike_total[wellIndex][erccIndex]++;
							if(!umi_seen.count(umicode)){
								umi_seen.insert(umicode);
								spike_umi[wellIndex][erccIndex]++;	
							}
						}
						if(!umi_seen_mm.count(umicode)){
							umi_seen_mm.insert(umicode);
							spike_umi_mm[wellIndex][erccIndex]++;	
						}						
					}
					else{
						const unsigned int geneIndex=category;
						total_mm[wellIndex][geneIndex]++;
						if (!multiGeneHit){
							total[wellIndex][geneIndex]++;
							if(!umi_seen.count(umicode)){
								umi_seen.insert(umicode);
								umi[wellIndex][geneIndex]++;	
							}
						}
						if(!umi_seen_mm.count(umicode)){
							umi_seen_mm.insert(umicode);
							umi_mm[wellIndex][geneIndex]++;	
						}
					}		
				}	
				fclose(fp);
			}
		}
	}
	
	bool checkPush(vector<string> items, MapPosition &best_hit_list){
		//check if it has been seen
		//>= is used to check MAX_BEST because the best hit is also in the list unlike the code for bwa where the extra hits are in the list
		//fprintf(stderr,"%d %d : ",best_hit_list.gene.size(),MAX_BEST);
		//fprintf(stderr,"%s %d %d : ",items[14].c_str(),stoi(splitStrIndex(items[14],":",-1)),MAX_EDIT_DISTANCE);
		//fprintf(stderr,"%s %d\n",items[9].c_str(),polyACheck(items[9]));
		
		if(best_hit_list.gene.size() >= MAX_BEST  || stoi(splitStrIndex(items[14],":",-1)) > MAX_EDIT_DISTANCE || polyACheck(items[9])) return 0;
		return best_hit_list.insert(items[2],stoi(items[3]));
	}	
	
	bool flushList(vector<string> items,string &barcode,unsigned int &barcodeIndex,MapPosition &best_hit_list){
		best_hit_list.clear(); 
		vector <string> tempItems;
		splitStr(items[0],":",tempItems);
		string fullBarcode=tempItems[tempItems.size()-1];
		barcode=fullBarcode.substr(6,UMISIZE); //change if barcode size changes
		barcodeIndex=0;
		if(ambigCheck(barcode,barcodeIndex)) return 0;
		return 1;
	} 
	~Counts(){
		delete[] spike_total[0];
		delete[] spike_umi[0];   
		delete[] total[0];
		delete[] umi[0];   
		delete[] spike_total_mm[0];
		delete[] spike_umi_mm[0];   
		delete[] total_mm[0];
		delete[] umi_mm[0];
		delete[] total_reads;
		delete[] assigned_reads;
		delete[] assigned_aligned_reads;	
		delete[] assigned_mito_reads;
		delete[] assigned_mito_umi;
		delete[] assigned_unknown_reads;
		delete[] assigned_unknown_umi;
		delete[] total_reads_mm;
		delete[] assigned_reads_mm;
		delete[] assigned_aligned_reads_mm;	
		delete[] assigned_mito_reads_mm;
		delete[] assigned_mito_umi_mm;
		delete[] assigned_unknown_reads_mm;
		delete[] assigned_unknown_umi_mm;
	}	
};	

using namespace std;   
int main(int argc, char *argv[]){
	int nbins=16;
	int binSize=0;
	string doneFile="";
	uint32_t posMask=(1 << nbins+1) -1;
	bool geneLevelFilter=0,filteredSAMfiles=0;
	string sample_id="",sym2ref="", ercc_fasta="", barcodes="", aligned_dir="", dge_dir="",countsFile="";
	int opt,verbose=0,nThreads=1;
	struct optparse options;
	optparse_init(&options, argv);	
 
	string errmsg="umimerge_parallel v?hfn:gd:o:t:i:s:e:m:c:b:a:p:\n-h -?  (display this message)\n-v (Verbose mode)\n-f Merge filtered SAM files (*.saf) instead of full SAM files\n-s <sym2ref file>\n-e <ercc_fasta file>\n-b <barcode_file>\n-a <aligns directory>\n-o <dge_dir(counts directory)>\n-c <binary file with umicounts - default is UMIcounts.bin in the input SAM directory>\n-t <number of threads (1)>\n-p <bin size for UMI position based filtering i.e 0 bits means reads with identical UMIs are discarded if they have same mapping position; 1 bit means reads with identical UMIs are discarded if their mapping position falls into same 2 basepair bin; 2 bit mean 4 basepair bins etc... \n\nRequired params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir\n\nExample:\n\numimerge_parallel -i RNAseq_20150409 -s  References/Broad_UMI/Human_RefSeq/refGene.hg19.sym2ref.dat -e References/Broad_UMI/ERCC92.fa -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -a Aligns -o Counts -t 4 -p 0\n";
	while ((opt = optparse(&options, "v?hfn:go:t:i:s:e:m:c:b:a:p:")) != -1) {
		switch (opt){
			case 'v':
				verbose=1;
			break;
			case 'd':
				doneFile=string(options.optarg);
			case 'n':
				nbins=atoi(options.optarg);
				posMask=(1 << nbins+1) -1;
			 //set nbins - default is 16 bits otherwise  
			break;
			case 'f':
			 //use filtered SAM files and not sam files
				filteredSAMfiles=1;
			break;
			case 'g':
			 //set gene level filtering - overrides other options
				geneLevelFilter=1;
			break;
			case 'o':
				dge_dir=string(options.optarg);
			break;
			case 't':
				nThreads=atoi(options.optarg);
			break;
			case 'i':
				sample_id=string(options.optarg);
			break;  
			case 's':
				sym2ref=string(options.optarg);
			break;
			case 'e':
				ercc_fasta=string(options.optarg);
			break;
			case 'c':
				countsFile=string(options.optarg);
			break;  
			case 'b':
				barcodes=string(options.optarg);
			break;
			case 'a':
				aligned_dir=string(options.optarg);
				if(countsFile ==""){
					countsFile=aligned_dir+"/"+"UMIcounts.bin";
				}	
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
	if(sample_id =="" || sym2ref=="" ||  ercc_fasta=="" || barcodes=="" || aligned_dir=="" || dge_dir==""){
		fprintf(stderr,"Required params are -i sample_id -s sym2ref -e ercc_fasta -b barcodes -a aligned_dir -o dge_dir\n");
		exit(EXIT_FAILURE);
	}	
	
	
	unordered_map<string,string> refseq_to_gene;
	vector<string>erccList, geneList,wellList;
	unordered_map<string,unsigned int>well_to_index;
	unordered_map<string,unsigned int>ercc_to_index;
	unordered_map<string,unsigned int>gene_to_index;
	vector<string> unknown_list;
	
	string plateID=readWells(barcodes,well_to_index,wellList);
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
		
	Counts count(nThreads,erccList.size(),geneList.size());
	//check if doneFile is to be used for asynchronous operation
	if doneFile != "") doneFile=aligned_dir+"/"+doneFile
	
	if (filteredSAMfiles) count.merge_filter(doneFile,nThreads,erccList,geneList,wellList,aligned_dir,barcodePanel,refseq_to_gene,well_to_index,ercc_to_index,gene_to_index, posMask,binSize,nbins,geneLevelFilter);
	else count.merge_parallel(doneFile,nThreads,erccList,geneList,wellList,aligned_dir,barcodePanel,refseq_to_gene,well_to_index,ercc_to_index,gene_to_index, posMask,binSize,nbins,geneLevelFilter);
		count.print(dge_dir,sample_id,nThreads,erccList,geneList,wellList);
	return 1;		
	 
}
bool doneFileExists(const std::string &doneFile){
	if (FILE *fp = fopen(name.c_str(), "r")) {
		fclose(fp);
		return 0;
	} 
	else {
		return 1;
	} 		
}
