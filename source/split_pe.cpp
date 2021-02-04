#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

extern "C" {
 #include "optparse.h"  
}
bool checkNames(char *line1,char *line2);
bool getLines(gzFile fp1, gzFile fp2, char *bufferR1, char *bufferR2, bool *mismatch);
bool closeFile(FILE *fp,gzFile gzfp, std::string filename,bool writeFile);
bool writeFilename(std::string filename,std::string suffix);
bool openFiles(FILE **ofp1,FILE **ofp2, gzFile *ofpgz1, gzFile *ofpgz2, std::string R1stem,std::string R2stem, std::string outputDir, int numFiles,bool compressFlag);
int writeLines(FILE *ofp,gzFile ofpgz,char *buf,unsigned int buflen,bool compressFlag);

std::string errmsg="split_pe h?vdt:zs:o:q:\n-h -? (display this message)\n-v (Verbose mode)\n-z Compress the output as gzip files\n-s The maximum size of the split file in KB\n):\n-t <number of threads(1)>\n-o <Output Directory>\n<R1file.1> <R2file.1>..<R1file.N> <R2file.N>\n\nExample:\numisplit -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz\n";
  
int main(int argc, char *argv[]){
    bool compressFlag=0,writeDoneFiles=0;
    char verbose=0;
    uint64_t maxSizeKB=0;  
    int nameMismatch=0;
    
    char *arg=0,*outputFileName=0;
    struct optparse options;
    int opt;
    optparse_init(&options, argv);
    int nThreads=1;
    int filter=0;
    std::string outputDir="";
    std::vector<std::string> inputFiles;
 
 //parse flags
    while ((opt = optparse(&options, "h?vdt:zs:o:")) != -1) {
        switch (opt){
            case 'v':
                verbose=1;
			break;
			case 'd':
				writeDoneFiles=1;
			break;
            case 'z':
                compressFlag=1;
            break;
            case 's':
                maxSizeKB=std::stoull(options.optarg);
            break;     
			case 'o':
                outputDir=std::string(options.optarg);
			break;
            case 't':
                nThreads=atoi(options.optarg);
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
	if(outputDir==""){
		fprintf(stderr,"must give an output directory\n");
		exit(EXIT_FAILURE);
	}	
    //parse file arguments
    //these will be R1 followed by R2 arguments
    while ((arg = optparse_arg(&options))){
        inputFiles.push_back(std::string(arg));
	}
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"Output directory is %s\n",outputDir.c_str());
		fprintf(stderr,"Number of threads %d\n",nThreads);		
		if(writeDoneFiles)fprintf(stderr,"writing out a file to indicate that we have finished with writing\n");
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
	//this assumes linux style directories 
	//make a subDir to store splits because the demux version of this makes subdirs for each well
	std::string splitDir=outputDir+"/splits/";
    fs::create_directory(fs::system_complete(splitDir));

    #pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
		const uint64_t maxSize=1024*maxSizeKB;
		const int tid=omp_get_thread_num();
		int numFiles=0;
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
        std::string R1stem,R2stem;
        R1=R1.stem();
        R2=R2.stem();
        while (R1.extension().string() == ".gz" || R1.extension().string() == ".fastq" || R1.extension().string() == ".fq")   
			R1=R1.stem();
		while (R2.extension().string() == ".gz" || R2.extension().string() == ".fastq" || R2.extension().string() == ".fq")   
			R2=R2.stem();   
        R1stem=R1.stem().string();
        R2stem=R2.stem().string();
        FILE *ofp1=0,*ofp2=0;
        gzFile ofpgz1=0,ofpgz2=0;
        uint64_t fileSize1,fileSize2;
        char linesR1[4096], linesR2[4096];
        bool mismatch=0;
        openFiles(&ofp1,&ofp2,&ofpgz1,&ofpgz2,R1stem,R2stem,splitDir,numFiles,compressFlag);
        while (getLines(fp1,fp2,linesR1,linesR2,&mismatch)){
            if (mismatch){
                nameMismatch++;
		 		if(verbose)fprintf(stderr,"Warning - mismatch of names in R1 and R2\n");
            }
            const unsigned int len1=strlen(linesR1);
            const unsigned int len2=strlen(linesR2);
			fileSize1+=len1;
			fileSize2+=len2;
            if(maxSize && (fileSize1 > maxSize || fileSize2 > maxSize)){
				openFiles(&ofp1,&ofp2,&ofpgz1,&ofpgz2,R1stem,R2stem,splitDir,numFiles,compressFlag);
				fileSize1=len1;
				fileSize2=len2;
				numFiles++;
			}
			writeLines(ofp1,ofpgz1,linesR1,len1,compressFlag);
			writeLines(ofp2,ofpgz2,linesR2,len2,compressFlag);
		}
		std::string outputFileR1=outputDir+"/"+R1stem;
		std::string outputFileR2=outputDir+"/"+R2stem;
		closeFile(0,fp1,outputFileR1,writeDoneFiles);
		closeFile(0,fp2,outputFileR2,writeDoneFiles);		
	}
    return 0;  
}
int writeLines(FILE *ofp,gzFile ofpgz,char *buf,unsigned int buflen,bool compressFlag){
	if (compressFlag){
	   return gzwrite(ofpgz,buf,buflen);
	}
    else return fwrite(buf,1,buflen,ofp);
}
bool getLines(gzFile fp1, gzFile fp2, char *bufferR1, char *bufferR2, bool *mismatch){
    const int lineSize=1024;
    char *bufptr1=bufferR1;
    char *bufptr2=bufferR2; 
    //check first line for mismatch
    if(!gzgets(fp1,bufferR1,lineSize) || !gzgets(fp2, bufferR2,lineSize)) return 0;
    *mismatch=checkNames(bufferR1,bufferR2);
    bufptr1+=strlen(bufptr1);
    bufptr2+=strlen(bufptr2);
    for (int i=0;i<3;i++){
        if(!gzgets(fp1,bufptr1,lineSize) || !gzgets(fp2, bufptr2,lineSize)) return 0;
		bufptr1+=strlen(bufptr1);
		bufptr2+=strlen(bufptr2);
	}
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

bool openFiles(FILE **ofp1,FILE **ofp2, gzFile *ofpgz1, gzFile *ofpgz2, std::string R1stem,std::string R2stem, std::string outputDir, int numFiles,bool compressFlag){
	const int numFieldSize=6;
    std::string numFilesString=std::to_string(numFiles);
    if (numFilesString.length() < numFieldSize) numFilesString = std::string(numFieldSize - numFilesString.length(), '0') + numFilesString;
	std::string outputfileR1 =outputDir+"/"+R1stem+"_"+ numFilesString +".fq";
	std::string outputfileR2 =outputDir+"/"+R2stem+"_"+ numFilesString +".fq";
	if (compressFlag){
	  outputfileR1+=".gz";
	  outputfileR2+=".gz";
	  if (*ofpgz1)gzclose(*ofpgz1);
	  if (*ofpgz2)gzclose(*ofpgz2);
	  *ofpgz1 = gzopen(outputfileR1.c_str(), "wb");
	  *ofpgz2 = gzopen(outputfileR2.c_str(), "wb");
	  if (!*ofpgz1 || !*ofpgz2) return 0;
	}
	else{
	  if (*ofp1)fclose(*ofp1);
	  if (*ofp2)fclose(*ofp2);
	  *ofp1 = fopen(outputfileR1.c_str(), "w");
	  *ofp2 = fopen(outputfileR2.c_str(), "w");
	  if (!*ofp1 || !*ofp2) return 0;		
    }
    return 1;
}
bool closeFile(FILE *fp,gzFile gzfp, std::string filename,bool writeFile){
	if (fp) fclose(fp);
	else if (gzfp) gzclose(gzfp);
	else return 1;
	if (writeFile) writeFilename(filename,"done");
	return 0;	
}
bool writeFilename(std::string filename,std::string suffix){
	std::string outputFilename=filename+"."+suffix;
	FILE *fp = fopen(outputFilename.c_str(),"w");
	if (fp){
		fprintf(fp,"%s",filename.c_str());
		fclose(fp);
		return 0;
	}
	return 1;
}
