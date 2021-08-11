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

std::string errmsg="Put informative error message here\n";
  
int main(int argc, char *argv[]){
    bool decompressFlag=0,writeDoneFiles=0;
    char verbose=0;
    uint64_t maxLines=2;
    std::string basename="part";
    std::string suffix="fq";
    char *arg=0,*outputFileName=0;
    struct optparse options;
    int opt;
    optparse_init(&options, argv);
    int nThreads=1;
    std::string outputDir="splitFiles";
    std::vector<std::string> inputFiles;
 
 //parse flags
    while ((opt = optparse(&options, "h?vdn:m:t:zs:o:")) != -1) {
        switch (opt){
            case 'v':
                verbose=1;
			break;
			case 'd':
				writeDoneFiles=1;
			break;
            case 'z':
                decompressFlag=1;
            break;
            case 'n':
                basename=std::string(options.optarg);
			break;
			case 'm':
                maxLines=atoi(options.optarg);
			break;     
			case 'o':
                outputDir=std::string(options.optarg);
			break;
			case 's':
                suffix=std::string(options.optarg);
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
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"Output directory is %s\n",outputDir.c_str());
		fprintf(stderr,"Number of threads %d\n",nThreads);		
		if(writeDoneFiles)fprintf(stderr,"writing out a file to indicate that we have finished with writing\n");
		int i=0;	
	}
	//this assumes linux style directories 
	std::string splitDir=outputDir+"/";
    fs::create_directory(fs::system_complete(splitDir));
//  #pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	char buffer[8192];
	char name[8192];
	char partName[8192];
	int filenumber=0;

	int nlines=0;
	if (decompressFlag){
		fprintf(stderr,"gunzipping\n",name);
		sprintf (name,"%s%s_%06d_.%s",splitDir.c_str(),basename.c_str(),filenumber,suffix.c_str());
		sprintf (partName,"%s%s_%06d_.part",splitDir.c_str(),basename.c_str(),filenumber);
		FILE *ofp = fopen(partName, "w");
		gzFile ifp = gzopen("/dev/stdin","r");
		while (gzgets(ifp,buffer,sizeof(buffer))){
			nlines++;
			if 	(nlines > maxLines){
				nlines=1;
				filenumber++;
				fprintf(stderr,"writing %s\n",name);
				std::rename(partName,name);
				sprintf (name,"%s%s_%06d_.%s",splitDir.c_str(),basename.c_str(),filenumber,suffix.c_str());
				sprintf (partName,"%s%s_%06d_.part",splitDir.c_str(),basename.c_str(),filenumber);
				ofp = fopen(partName, "w");
			}		  
			fputs(buffer,ofp);
		}
		if (ofp) fclose(ofp);
		fprintf(stderr,"writing %s\n",name);
		std::rename(partName,name);
	}		
    else{
		sprintf (name,"%s%s_%06d_.%s",splitDir.c_str(),basename.c_str(),filenumber,suffix.c_str());
		fprintf(stderr,"writing %s\n",name);
		FILE *ofp = fopen(name, "w");
		while (fgets(buffer,sizeof(buffer),stdin)){
			nlines++;
			if 	(nlines > maxLines){
				nlines=1;
				filenumber++;
				fprintf(stderr,"writing %s\n",name);
				std::rename(partName,name);
				sprintf (name,"%s%s_%06d_.%s",splitDir.c_str(),basename.c_str(),filenumber,suffix.c_str());
				sprintf (partName,"%s%s_%06d_.part",splitDir.c_str(),basename.c_str(),filenumber);
				ofp = fopen(partName, "w");
			}		  
			fputs(buffer,ofp);
		}
		if (ofp) fclose(ofp);
		fprintf(stderr,"writing %s\n",name);
		std::rename(partName,name);
	}
	if (writeDoneFiles){
		std::string doneName=splitDir+basename+"_done";
		FILE *ofp = fopen(doneName.c_str(), "w");
		fclose(ofp);
	}
	return 0;
}

