#!/bin/bash

inpFolder=/data/input
outFolder=/data/output

for i in `ls $inpFolder/*`; do
	echo $i
	a=$(basename $i);
	b=`echo $a | awk -F '.' '{print $1}'`;
	c=`echo $b | awk -F '-' '{print $1"-"$2"-"$3}'`;
	echo $a, $b, $c
	mainSegmentFeatures -t wsi -i $i -p $outFolder/$b -d 4096
	mainAggregateFeatures.py -i "$outFolder/$b*-features.csv" -p $c -o $outFolder/$c-aggregate.json -t json "$@" 
done
