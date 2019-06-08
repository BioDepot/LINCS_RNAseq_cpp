#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include "umitools.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

extern "C" {
 #include "optparse.h"  
}
bool checkNames(char *line1,char *line2);
bool Ncheck(const char *seq, const int size);
bool getLines(gzFile fp1, gzFile fp2, char *bufferR1[], char *bufferR2, bool *mismatch);
bool closeFile(FILE *fp,gzFile gzfp, string filename,bool writeFile);
bool writeFilename(string filename,string suffix);


#define R1_LENGTH 16 //default UMI length

string errmsg="umisplit h?vftd:m:N:o:b:l:q:\n-h -? (display this message)\n-v (Verbose mode)\n-f (filter and discard reads with ambiguous UMI and barcode (default is to keep))\n-z Compress the output as gzip files\n-s The maximum size of the split file in KB\n\n-l <Length of UMI (16):\n-t <number of threads(1)>\n-q <minimum quality of UMI base pair before changed to N (10)>-b <barcode file>\n-m <maximum number of mismatches tolerated in barcode (0)>\n-c <the file to output the umicounts to - default is UMIcounts.bin in the output directory>\n-N <maximum number of ambiguous base pairs tolerated in barcode (0)>\n-o <Output Directory>\n<R1file.1> <R2file.1>..<R1file.N> <R2file.N>\n\nExample:\numisplit -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz\n";
  
int main(int argc, char *argv[]){
    bool compressFlag=0,writeDoneFiles=0;
    char verbose=0,barcode=0;
    uint64_t maxSizeKB=0;  
    int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,UMILength=R1_LENGTH,nArgs=0;
    unsigned int barcodeSize=0;
    int adjustedQ=minQual+33;
    char *arg=0,*outputFileName=0;
    string barcodeFileName,countsFile="";
    struct optparse options;
    int opt;
    int mismatchTol=0,NTol=0;
    optparse_init(&options, argv);
    int nThreads=1;
    int filter=0;
    string outputDir="";
    vector<string> inputFiles;
 
 //parse flags
    while ((opt = optparse(&options, "h?vdft:zs:m:N:o:b:l:q:c:")) != -1) {
        switch (opt){
            case 'v':
                verbose=1;
			break;
			case 'd':
				writeDoneFiles=1;
			break;
			case 'f':
                filter=1;
            break;
            case 'z':
                compressFlag=1;
            break;
            case 'c':
                countsFile=string(options.optarg);
            break; 
            case 's':
                maxSizeKB=std::stoull(options.optarg);
            break;     
			case 'o':
                outputDir=string(options.optarg);
			break;
            case 'm':
                mismatchTol=atoi(options.optarg);
            break;
            case 't':
                nThreads=atoi(options.optarg);
            break;         
            case 'N':
                NTol=atoi(options.optarg);
            break;
            case 'b':
                barcodeFileName=options.optarg;
            break;
            case 'l':
                UMILength=atoi(options.optarg);//not the well barcode
            break;
            case 'q':
                minQual=atoi(options.optarg);
                adjustedQ=33+minQual;
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
	if(outputDir=="" || barcodeFileName==""){
		fprintf(stderr,"must give an output directory and barcodeFilename\n");
		exit(EXIT_FAILURE);
	}	
    //parse file arguments
    //these will be R1 followed by R2 arguments
    while ((arg = optparse_arg(&options))){
        inputFiles.push_back(string(arg));
	}
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"UMI Length %d\n",UMILength);
		fprintf(stderr,"Barcode fileName %s\n",barcodeFileName.c_str());
		fprintf(stderr,"Output directory is %s\n",outputDir.c_str());
		fprintf(stderr,"Minimum quality %d\n",minQual);		
		fprintf(stderr,"Number of threads %d\n",nThreads);		
		fprintf(stderr,"Maximum number of mismatches in barcode %d\n",mismatchTol);		
		fprintf(stderr,"Maximum number of unknown bases in barcode %d\n",NTol);
		if(filter)fprintf(stderr,"filtering out unmatched barcodes \n");
		if(writeDoneFiles)fprintf(stderr,"writing out a file to indicate that we have finished with writing\n");
    else (stderr,"saving unmatched barcodes to 'X' directory \n");
		int i=0;
		while (i<inputFiles.size()){
			fprintf(stderr,"R1 file %s\n",inputFiles[i++].c_str());
			if(i == inputFiles.size()){
				fprintf(stderr,"missing corresponding R2 file\n");
				exit(EXIT_FAILURE);
			}
			fprintf(stderr,"R2 file %s\n",inputFiles[i++].c_str());		
		}		
	}
    fs::create_directory(fs::system_complete(outputDir));
    
    #if NWELLS > 256
        umipanel<uint32_t,uint16_t> **barcodePanel=new umipanel<uint32_t,uint16_tr>*[nThreads];
        for(int i=0;i<nThreads;i++)
            barcodePanel[i]=new umipanel<uint32_t,uint16_t> (barcodeFileName,mismatchTol,NTol);
	#else
        umipanel<uint32_t,unsigned char> **barcodePanel=new umipanel<uint32_t,unsigned char>*[nThreads];
        for(int i=0;i<nThreads;i++)
            barcodePanel[i]=new umipanel<uint32_t,unsigned char> (barcodeFileName,mismatchTol,NTol);
	#endif
	
	//create directories for each well
	for(int i=0;i<NWELLS;i++){
		auto p=fs::path(outputDir+"/"+barcodePanel[0]->wells[i]);
		fs::create_directory(fs::system_complete(p));
	}	
	//create bad directory	
	if(!filter){
        auto p=fs::path(outputDir+"/X");
        fs::create_directory(fs::system_complete(p));
	}
    #pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
		const uint64_t maxSize=1024*maxSizeKB;
		const int tid=omp_get_thread_num();
		const int wellSequenceSize=barcodePanel[tid]->barcodeSize;
		gzFile fp1=0, fp2=0;  
		fp1 = gzopen(inputFiles[i].c_str(), "r");
		if(fp1 == Z_NULL){
			fprintf(stderr,"unable to open %s\n",inputFiles[i].c_str());
			exit(EXIT_FAILURE);
		}
		fp2 = gzopen(inputFiles[i+1].c_str(), "r");	
		if(fp2 == Z_NULL){
			fprintf(stderr,"unable to open %s\n",inputFiles[i+1].c_str());
			exit(EXIT_FAILURE);
		}
  
        fs::path R1(inputFiles[i]);
        fs::path R2(inputFiles[i+1]);
        string R1stem,R2stem;
        R1=R1.stem();
        R2=R2.stem();
        while (R1.extension().string() == ".gz" || R1.extension().string() == ".fastq" || R1.extension().string() == ".fq")   
			R1=R1.stem();
		while (R2.extension().string() == ".gz" || R2.extension().string() == ".fastq" || R2.extension().string() == ".fq")   
			R2=R2.stem();   
        R1stem=R1.stem().string();
        R2stem=R2.stem().string();
        FILE *ofp,*ofps[NWELLS+1];
        vector<string> filenames(NWELLS+1);
        gzFile ofpgz,ofpsgz[NWELLS+1];
        uint64_t fileSizes[NWELLS+1];
        uint16_t numberofFiles[NWELLS+1];
        memset(fileSizes,0,sizeof(fileSizes));
        memset(numberofFiles,0,sizeof(numberofFiles));
  
        for(int j=0;j<NWELLS+1;j++){
			//zeroing is necessary - we use the zero value in the file pointer instead of passing the compress flag to know which type to use for opening/closing
			ofps[j]=0;
			ofpsgz[j]=0;
			string file;
			if(!filter && !j){
				file=outputDir+"/X"+"/"+R1stem+"R2_"+"X_"+std::to_string(numberofFiles[0])+".fq";
				numberofFiles[0]++;
				if (compressFlag){
					file=file+".gz";
					ofpsgz[0]=gzopen(file.c_str(),"wb");			
				}
				else{
                    ofps[0]=fopen(file.c_str(),"w");			
                    if(!ofps[j]){
                        fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
                    }
				}
				filenames[j]=file;	
			}	
			else if(j){
                file =outputDir+"/"+barcodePanel[tid]->wells[j-1]+"/"+R1stem+"R2_"+barcodePanel[tid]->wells[j-1]+"_"+std::to_string(numberofFiles[j])+".fq";
				if (compressFlag){
					file=file+".gz";
					ofpsgz[j]=gzopen(file.c_str(),"wb");			
				}
				else{
                    ofps[j]=fopen(file.c_str(),"w");
                    if(!ofps[j]){
                        fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
                    }	
                }
                filenames[j]=file;	
                numberofFiles[j]++;
            }
        }
        char linesR1buf[4096];
        char *linesR1[4]={linesR1buf,linesR1buf+1024,linesR1buf+2048,linesR1buf+3072},linesR2[1024];
        bool mismatch=0;
        while (getLines(fp1,fp2,linesR1,linesR2,&mismatch)){
            if (mismatch){
                nameMismatch++;
		 		if(verbose)fprintf(stderr,"Warning - mismatch of names in R1 and R2\n");
            }
            //adjust quality and find wells
            int k=0;
            char *cptr=linesR1[1];
            char *qptr=linesR1[3];
			while(*cptr && k < UMILength){
                if(*qptr<adjustedQ)*cptr='N';
                cptr++;qptr++;k++;
            } 
            //check if there is a N
            if(filter && Ncheck(linesR1[1]+wellSequenceSize,UMILength-wellSequenceSize)) continue;
            const unsigned int barcodeIndex=barcodePanel[tid]->bestMatch(linesR1[1]);
            if(filter && !barcodeIndex)continue;
			//skip if ambiguous barcode and get unique barcode index from sequence
			string umiBarcode(linesR1[1]+wellSequenceSize,UMILength-wellSequenceSize);
			unsigned int  umiBarcodeIndex=0;
			if (filter && ambigCheck(umiBarcode,umiBarcodeIndex)) continue;
			ofp=0;ofpgz=0;
			if(compressFlag) ofpgz=ofpsgz[barcodeIndex];
			else ofp=ofps[barcodeIndex];
            //check maxSize when writing
            if(maxSize){
                uint16_t lineBytes=strlen(linesR2)+strlen(linesR1[0])+UMILength+7;
                if (barcodeIndex && lineBytes+fileSizes[barcodeIndex] > maxSize){
					closeFile(ofp,ofpgz,filenames[barcodeIndex],writeDoneFiles);
					string file =outputDir+"/"+barcodePanel[tid]->wells[barcodeIndex-1]+"/"+R1stem+"R2_"+barcodePanel[tid]->wells[barcodeIndex-1]+"_"+std::to_string(numberofFiles[barcodeIndex])+".fq";
					if (compressFlag){
					 file=file+".gz";
					 ofpsgz[barcodeIndex]=gzopen(file.c_str(),"wb");	
					 ofpgz=ofpsgz[barcodeIndex];
					}
					else{
                        ofps[barcodeIndex]=fopen(file.c_str(),"w");
                        ofp=ofps[barcodeIndex];
					}
                    numberofFiles[barcodeIndex]++;
                    filenames[barcodeIndex]=file;
					fileSizes[barcodeIndex]=0;
				}
                fileSizes[barcodeIndex]+=lineBytes;
			}
            {
				char buffer[1024];
				memset(buffer,0,sizeof(buffer));
				char *bufPtr=buffer;
				char *dest=linesR1[0];
                while(*dest){
					*bufPtr++=*dest++;	
				}
                *bufPtr++=':';
				memcpy(bufPtr,linesR1[1],UMILength);
				bufPtr+=UMILength;
				*bufPtr++='\n';
				dest=linesR2;
				while(*dest){
					*bufPtr++=*dest++;	
				}
                if (compressFlag) gzwrite(ofpgz,buffer,bufPtr-buffer);
                else fwrite(buffer,1,bufPtr-buffer,ofp);
			}
		}
		//if not filtered we write the unreadable barcodes to another file
		if(!filter)closeFile(ofps[0],ofpsgz[0],"",0);
		 
		for(int j=1;j<NWELLS+1;j++){
			closeFile(ofps[j],ofpsgz[j],filenames[j],writeDoneFiles);
		}
		string outputFileR1=outputDir+"/"+R1stem;
		string outputFileR2=outputDir+"/"+R2stem;
		closeFile(0,fp1,outputFileR1,writeDoneFiles);
		closeFile(0,fp2,outputFileR2,writeDoneFiles);		
	}
    delete[] barcodePanel;
    return 0;  
}
bool getLines(gzFile fp1, gzFile fp2, char *bufferR1[], char *bufferR2, bool *mismatch){
    const int lineSize=1024;
    if(!gzgets(fp1,bufferR1[0],lineSize) || !gzgets(fp2, bufferR2,lineSize)) return 0;
    //checkNames also terminates bufferR1[0] at first space
    *mismatch=checkNames(bufferR1[0],bufferR2);
    if (*mismatch)fprintf(stderr,"warning mismatch %s %s",bufferR1[0],bufferR2);
        
    if (!gzgets(fp1, bufferR1[1], lineSize)) return 0;
    if (!gzgets(fp1, bufferR1[2], lineSize)) return 0; 
    if (!gzgets(fp1, bufferR1[3], lineSize)) return 0;
    char *bufptr=bufferR2;
    if (!gzgets(fp2, bufferR2,lineSize)) return 0;
    bufptr=bufferR2+strlen(bufferR2);
    if (!gzgets(fp2, bufptr,lineSize)) return 0;
    bufptr+=strlen(bufptr);
    if (!gzgets(fp2, bufptr, lineSize)) return 0;
    return 1;    
}
bool checkNames(char *line1,char *line2){
    char *a=line1,*b=line2;
    while(*a && *b){
        if (*a ==' ' && *b== ' '){
            *a=0;
            return 0;
        }        
        if (*a != *b) return 1;
        a++;
        b++;
    }
    return 0;    
}

bool Ncheck(const char *seq, const int size){
	for(int i=0;i<size;i++)
	 if(seq[i] == 'N')return 1;
	return 0;
}
bool closeFile(FILE *fp,gzFile gzfp, string filename,bool writeFile){
	if (fp) fclose(fp);
	else if (gzfp) gzclose(gzfp);
	else return 1;
	if (writeFile) writeFilename(filename,"done");
	return 0;	
}
bool writeFilename(string filename,string suffix){
	string outputFilename=filename+"."+suffix;
	FILE *fp = fopen(outputFilename.c_str(),"w");
	if (fp){
		fprintf(fp,"%s",filename.c_str());
		fclose(fp);
		return 0;
	}
	return 1;
}
