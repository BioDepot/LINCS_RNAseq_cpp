#!/bin/bash

sudo docker build -t rnaseq-umi-cpp -f Dockerfile.build  ${PWD}
sudo docker run --rm -v ${PWD}:/local rnaseq-umi-cpp /bin/sh -c "cp -r source/w* /local/."
sudo docker run --rm -v ${PWD}:/local biodepot/alpine-bwa:3.9.2-0.7.15 /bin/sh -c "cp -r /usr/local/bin/bwa /local/."
sudo docker build -t biodepot/rnaseq-umi-cpp:profiler  -f Dockerfile-profiler  ${PWD}
sudo rm -rf w96 w384
