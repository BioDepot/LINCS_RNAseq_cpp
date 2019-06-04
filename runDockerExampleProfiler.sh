#!/bin/bash
./download.sh --decompress --directory ./data https://drive.google.com/open?id=16QHgiI_9QYuCukjZQmw03u7w374R3kt1 
echo "./scripts/docker_fast_run-alignment-analysis_profiler.sh ${PWD}"
./scripts/docker_fast_run-alignment-analysis_profiler.sh ${PWD}
