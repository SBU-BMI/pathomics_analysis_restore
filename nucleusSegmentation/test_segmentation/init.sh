#!/bin/bash

# Build image from dockerfile
docker build -t sbubmi/test_segmentation .

# Start container from image
docker_image="sbubmi/test_segmentation"
docker_name="$USER-test_segmentation"
python ../script/run_docker_segment.py start $docker_name $docker_image
