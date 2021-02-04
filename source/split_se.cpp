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
bool getLines(gzFile fp1, char *bufferR1);
bool closeFile(FILE *fp,gzFile gzfp, std::string filename,bool writeFile);
bool writeFilename(std::string filename,std::string suffix);
bool openFiles(FILE **ofp1, gzFile *ofpgz1, std::string R1stem,std::string outputDir, int numFiles,bool compressFlag);
int writeLines(FILE *ofp,gzFile ofpgz,char *buf,unsigned int buflen,bool compressFlag);

std::string errmsg="split_pe h?vdt:zs:o:q:\n-h -? (display this message)\n-v (Verbose mode)\n-z Compress the output as gzip files\n-s The maximum size of the split file in KB\n):\n-t <number of threads(1)>\n-o <Output Directory>\n<R1file.1> <R2file.1>..<R1file.N> <R2file.N>\n\nExample:\numisplit -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz\n";
  
int main(int argc, char *argv[]){
    bool compressFlag=0,writeDoneFiles=0;
    char verbose=0;
    uint64_t maxSizeKB=0;  
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
			fprintf(stderr,"fastq file %s\n",inputFiles[i++].c_str());	
		}		
	}
	//this assumes linux style directories 
	//make a subDir to store splits because the demux version of this makes subdirs for each well
	std::string splitDir=outputDir+"/splits/";
    fs::create_directory(fs::system_complete(splitDir));

    #pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i++){
		const uint64_t maxSize=1024*maxSizeKB;
		const int tid=omp_get_thread_num();
		int numFiles=0;
		gzFile fp1=0;  
		fp1 = gzopen(inputFiles[i].c_str(), "r");
		if(fp1 == Z_NULL){
			fprintf(stderr,"unable to open %s\n",inputFiles[i].c_str());
			exit(EXIT_FAILURE);
		}
        fs::path R1(inputFiles[i]);
        std::string R1stem;
        R1=R1.stem();
        while (R1.extension().string() == ".gz" || R1.extension().string() == ".fastq" || R1.extension().string() == ".fq")   
			R1=R1.stem();;   
        R1stem=R1.stem().string();
        FILE *ofp1=0;
        gzFile ofpgz1=0;
        uint64_t fileSize1;
        char linesR1[4096];
        openFiles(&ofp1,&ofpgz1,R1stem,splitDir,numFiles,compressFlag);
        while (getLines(fp1,linesR1)){
            const unsigned int len1=strlen(linesR1);
			fileSize1+=len1;
            if(maxSize && fileSize1 > maxSize){
				openFiles(&ofp1,&ofpgz1,R1stem,splitDir,numFiles,compressFlag);
				fileSize1=len1;
				numFiles++;
			}
			writeLines(ofp1,ofpgz1,linesR1,len1,compressFlag);
		}
		std::string outputFileR1=outputDir+"/"+R1stem;
		closeFile(0,fp1,outputFileR1,writeDoneFiles);
	}
    return 0;  
}
int writeLines(FILE *ofp,gzFile ofpgz,char *buf,unsigned int buflen,bool compressFlag){
	if (compressFlag){
	   return gzwrite(ofpgz,buf,buflen);
	}
    else return fwrite(buf,1,buflen,ofp);
}
bool getLines(gzFile fp1, char *bufferR1){
    const int lineSize=1024;
    char *bufptr1=bufferR1;
    //check first line for mismatch
    for (int i=0;i<4;i++){
        if(!gzgets(fp1,bufptr1,lineSize)) return 0;
		bufptr1+=strlen(bufptr1);
	}
    return 1;    
}

bool openFiles(FILE **ofp1,gzFile *ofpgz1, std::string R1stem,std::string outputDir, int numFiles,bool compressFlag){
	const int numFieldSize=6;
    std::string numFilesString=std::to_string(numFiles);
    if (numFilesString.length() < numFieldSize) numFilesString = std::string(numFieldSize - numFilesString.length(), '0') + numFilesString;
	std::string outputfileR1 =outputDir+"/"+R1stem+"_"+ numFilesString +".fq";
	if (compressFlag){
	  outputfileR1+=".gz";
	  if (*ofpgz1)gzclose(*ofpgz1);
	  *ofpgz1 = gzopen(outputfileR1.c_str(), "wb");
	  if (!*ofpgz1) return 0;
	}
	else{
	  if (*ofp1)fclose(*ofp1);
	  *ofp1 = fopen(outputfileR1.c_str(), "w");
	  if (!*ofp1) return 0;		
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
