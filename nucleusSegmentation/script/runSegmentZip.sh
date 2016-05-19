#!/bin/bash

echo "$#"

if [ "$#" -ne 16 ]; then
	echo "Usage: runSegmentTileZip.sh -i <input image> -o <output prefix> -a <analysis id> -c <cancer type> -l <tileLeftX,tileLeftY> -s <tileSizeX,tileSizeY> -p <patchSizeX,patchSizeY> -z <zipFilename>"
	exit 1;
fi

while [[ $# > 1 ]]
do
	key="$1"

	case $key in
		-i)
			inpFile="$2"
			shift # past argument
		;;
		-o)
			outPrefix="$2"
			shift # past argument
		;;
		-a)
			analysisId="$2"
			shift # past argument
		;;
		-c)
			cancerType="$2"
			shift # past argument
		;;
		-l)
			topLeft="$2"
			shift # past argument
		;;
		-s)
			tileSize="$2"
			shift # past argument
		;;
		-p)
			patchSize="$2"
			shift # past argument
		;;
		-z)
			zipFilename="$2"
			shift # past argument
		;;
		--default)
			DEFAULT=YES
		;;
		*)
			echo "Unknown option";
			exit 1;	# unknown option
		;;
	esac
	shift # past argument or value
done

fileName=$(basename $inpFile);
patientId=`echo $fileName | awk -F '-' '{print $1"-"$2"-"$3}'`;

mainSegmentFeatures -t onetile -i $inpFile -p $outPrefix -s $topLeft -b $tileSize -d $patchSize -a $analysisId -e $analysisId -v mask:img:overlay -z $zipFilename
mainAggregateFeatures.py -i "$outPrefix*-features.csv" -p $patientId -a $analysisId -o $patientId\-aggr.json -t json -s $cancerType -l TP 
zip -ujr $zipFilename $patientId\-aggr.json


