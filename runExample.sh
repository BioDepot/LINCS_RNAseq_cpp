#!/bin/bash
download.sh --decompress --directory ./data https://drive.google.com/open?id=16QHgiI_9QYuCukjZQmw03u7w374R3kt1 
echo "./scripts/fast_run-alignment-analysis.sh ${PWD}"
./scripts/fast_run-alignment-analysis.sh ${PWD}
