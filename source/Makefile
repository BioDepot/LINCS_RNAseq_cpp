CC=gcc
CPP=g++
LD=g++ -L /usr/local/lib
ifndef ARM64
 m64FLAG=-m64 
 sse3FLAG=-msse3
 sse4FLAG=-msse4
endif
ifdef STATIC 
 STATIC_FLAG = -static
endif
ifdef DEBUG
CFLAGS= -I./ -I./optparse -ggdb -I -Wall -fopenmp -std=c++11 ${STATIC_FLAG}
else 
 ifdef AVX
  CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp ${m64FLAG} -mavx  ${sse4FLAG} -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES) 
 else
  ifdef SSE2
   ifdef 32Bit
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   else
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp ${m64FLAG} -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   endif
  else
   ifdef 32Bit
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp ${sse3FLAG} -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   else
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp ${m64FLAG} ${sse3FLAG} -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   endif 
  endif
 endif
endif
LIBS=-lboost_filesystem -lboost_system -lgomp -lm -lc -lgcc -lrt -lz
MERGELIBS= -lm -lc -lgcc -lrt -lz
MKDIR_P = mkdir -p
#ifdef ALPINE
 ALPINE_FLAG= -fPIC -static
#endif
OCFLAGS = -std=c99 -Wall -Wextra -g3 ${ALPINE_FLAG}
LFLAGS = -static-libstdc++

.PHONY: clean cleanDeps cleanBins all all96 all384 

default: all
all: all96 all384 clean
all96: clean directories96 umimerge_parallel umimerge_filter umisplit split_pe split_se streamSplit
all384: clean directories384 umimerge_parallel umimerge_filter umisplit split_pe split_se streamSplit

OBJ_umisplit = umisplit.o optparse.o 
OBJ_split_pe = split_pe.o optparse.o 
OBJ_streamSplit = streamSplit.o optparse.o 
OBJ_split_se = split_se.o optparse.o 
OBJ_umimerge_parallel = umimerge_parallel.o optparse.o 
OBJ_umimerge_filter = umimerge_filter.o optparse.o 
OBJ_merge_filter = merge_filter.o optparse.o 

directories96: 
	$(eval DIR+=w96/) 
	${MKDIR_P} ${DIR} 
directories384: 
	$(eval DIR+=w384/) 
	$(eval NWELLS+=w384/) 
	${MKDIR_P} ${DIR} 

umimerge_parallel : $(OBJ_umimerge_parallel)
	$(LD) -o ${DIR}umimerge_parallel $(LFLAGS) $(OBJ_umimerge_parallel) $(LIBS)
umimerge_parallel.o : umimerge_parallel.cpp 
	$(CPP) $(CFLAGS) -c umimerge_parallel.cpp			
	
umimerge_filter : $(OBJ_umimerge_filter)
	$(LD) -o ${DIR}umimerge_filter -static-libstdc++ $(OBJ_umimerge_filter) $(MERGELIBS)
umimerge_filter.o : umimerge_filter.cpp 
	$(CPP) $(CFLAGS) -c umimerge_filter.cpp			

umisplit : $(OBJ_umisplit)
	$(LD) -o ${DIR}umisplit $(LFLAGS) $(OBJ_umisplit) $(LIBS)
umisplit.o : umisplit.cpp 
	$(CPP) $(CFLAGS) -c umisplit.cpp
	
streamSplit : $(OBJ_streamSplit)
	$(LD) -o ${DIR}streamSplit $(LFLAGS) $(OBJ_streamSplit) $(LIBS)
streamSplit.o : streamSplit.cpp 
	$(CPP) $(CFLAGS) -c streamSplit.cpp

split_pe : $(OBJ_split_pe)
	$(LD) -o ${DIR}split_pe $(LFLAGS) $(OBJ_split_pe) $(LIBS)
split_pe.o : split_pe.cpp 
	$(CPP) $(CFLAGS) -c split_pe.cpp
	
split_se : $(OBJ_split_se)
	$(LD) -o ${DIR}split_se $(LFLAGS) $(OBJ_split_se) $(LIBS)
split_se.o : split_se.cpp 
	$(CPP) $(CFLAGS) -c split_se.cpp
								
optparse.o : optparse/optparse.c
	$(CC) $(OCFLAGS) -c optparse/optparse.c

clean: cleanBins cleanDeps

cleanDeps:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf optparse/*.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true

cleanBins: 
	-@rm  -f umimerge_parallel umisplit umimerge_filter split_pe split_se
