#!/bin/bash

#docker_image="sbubmi/test_segmentation"
#docker_name="$USER-test_segmentation"
docker_image="sbubmi/quip-jobs"
docker_name="quip-jobs"

# Build image from dockerfile
docker build -t $docker_image .

# Start container from image
python ../script/run_docker_segment.py start $docker_name $docker_image
