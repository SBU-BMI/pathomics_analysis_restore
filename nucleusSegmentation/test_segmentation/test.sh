#!/bin/bash

docker_name="$USER-test_segmentation"
cd ../script

# SEGMENT IMAGE
exec_id="20170112055937"
input_file="../test_segmentation/TCGA-CS-4938-01Z-00-DX1_12560_47520_500_500_LGG.png"

# One-liner:
#./run_docker_segment.py segment $docker_name $input_file $HOME/test_out.zip -t img -s 12560,47520 -b 500,500 -d 500,500 -a $exec_id -c TCGA-CS-4938-01Z-00-DX1 -p TCGA-CS-4938

# Line-breaks, easier to read:
./run_docker_segment.py \
segment \
$docker_name \
$input_file \
$HOME/test_out.zip \
-t img \
-j Y \
-s 12560,47520 \
-b 500,500 \
-d 500,500 \
-a $exec_id \
-c TCGA-CS-4938-01Z-00-DX1 \
-p TCGA-CS-4938
