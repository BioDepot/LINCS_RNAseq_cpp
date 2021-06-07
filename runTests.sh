#!/bin/bash
clean() {
 (cd $1; rm -rf Aligns/* Counts*/*)	
}

[ -d ./data ] || ./download.sh --decompress --directory ./data https://drive.google.com/open?id=16QHgiI_9QYuCukjZQmw03u7w374R3kt1 
[ -d ./source/w96 ] || (cd ./source; make)

#(cd ./source; rm -r w96; make clean; make)
clean ./data/LINCS
echo "creating sam files"
echo "./testScripts/generateSAM.sh ${PWD}"
./testScripts/generateSAM.sh ${PWD}
echo "creating saf files"
echo "./testScripts/generateSAF.sh ${PWD}"
./testScripts/generateSAF.sh ${PWD}
