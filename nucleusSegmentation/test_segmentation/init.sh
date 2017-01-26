#!/bin/bash

docker build -t sbubmi/test_segmentation .

git clone https://github.com/SBU-BMI/pathomics_analysis.git

cd pathomics_analysis/nucleusSegmentation/script

docker_image="sbubmi/test_segmentation"
docker_name="$USER-test_segmentation"
./run_docker_segment.py start $docker_name $docker_image
