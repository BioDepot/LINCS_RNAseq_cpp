#!/bin/bash
clean() {
 (cd $1; rm -rf Aligns/* Counts*/*)	
}

[ -d ./data ] || ./download.sh --decompress --directory ./data https://drive.google.com/open?id=16QHgiI_9QYuCukjZQmw03u7w374R3kt1 
[ -d ./source/w96 ] || (cd ./source; make)

#(cd ./source; rm -r w96; make clean; make)
clean ./data/LINCS
echo "testing on pe files"
echo "./testScripts/generatePE.sh ${PWD}"
./testScripts/generatePE.sh ${PWD} 
#clean ./data/LINCS
#echo "testing on non-multiplexed files"
#echo "./testScripts/generateNoMultiplex.sh ${PWD}"
#./testScripts/generateNoMultiplex.sh ${PWD} 
#clean ./data/LINCS
#echo "creating no-umi files"
#cho "./testScripts/generateNoUMI.sh ${PWD}"
#./testScripts/generateNoUMI.sh ${PWD} ./data/LINCS
#clean ./data/LINCS
#echo "creating sam_saf files"
#echo "./testScripts/generateSAMSAF.sh ${PWD}"
#./testScripts/generateSAMSAF.sh ${PWD}
#echo "creating sam files"
#echo "./testScripts/generateSAM.sh ${PWD}"
#./testScripts/generateSAM.sh ${PWD}
#echo "creating saf files"
#echo "./testScripts/generateSAF.sh ${PWD}"
#./testScripts/generateSAF.sh ${PWD}
#echo "creating old new sam files"
#echo "./testScripts/generateOldNewSAM.sh ${PWD}"
#./testScripts/generateOldNewSAM.sh ${PWD}
