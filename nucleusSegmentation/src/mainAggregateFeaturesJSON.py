#! /usr/bin/env python

import os.path
import glob
import sys, getopt
import socket
import numpy as np

def readCSV(dataDir):
    filePathnameList = glob.glob(dataDir) 
    allContent = [] 
    numFeatures = 0
    for j in filePathnameList: 
       print j
       f = open(j,'r')
       readContent = f.read().split('\n')
       numFeatures = len(readContent[0].split(','))
       readContent.pop(0)
       allContent += readContent 
       mydata = allContent.pop()
       f.close()   
    return (allContent,numFeatures)

def computeAggregateFeatures(allContent,numberOfFeatures):
    nObjects = len(allContent)
    allFeatureMatrix = np.empty([nObjects, numberOfFeatures], dtype=float)

    for iObject in range(nObjects):
        featuresOfThisObject = allContent[iObject]
        featuresOfThisObject = featuresOfThisObject.split(',')
        featuresOfThisObject = featuresOfThisObject[:numberOfFeatures]
        featuresOfThisObject = [float(x) for x in featuresOfThisObject]

        allFeatureMatrix[iObject, :] = np.asarray(featuresOfThisObject)

    [featureQ25, featureQ75, featureMedian] = np.percentile(allFeatureMatrix, [25, 50, 75], axis=0)
    subjectFeatures = np.concatenate((featureQ25, featureQ75, featureMedian))
    return subjectFeatures

def main(argv):
    dataDir = ''
    patientId = ''
    analysisId = ''
    outFile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:p:a:o:")
    except getopt.GetoptError:
       print 'mainAggregateFeatures.py [-h] [-i <inputfile pattern> -p <patient id> -a <analysis id> -o <output file>]'
       sys.exit(2)
    if len(argv) == 0:
       print 'mainAggregateFeatures.py [-h] [-i <inputfile pattern> -p <patient id> -a <analysis id> -o <output file>]'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'mainAggregateFeatures.py [-h] [-i <inputfile pattern> -p <patient id> -a <analysis id> -o <output file>]'
          sys.exit()
       elif opt in ("-i", "--ifile"):
          dataDir = arg
       elif opt == '-p':
          patientId = arg
       elif opt == '-a':
          analysisId = arg
       elif opt == '-o':
          outFile = arg

    if dataDir=='' or patientId=='' or analysisId=='' or outFile=='':
       print 'mainAggregateFeatures.py [-h] [-i <inputfile pattern> -p <patient id> -a <analysis id> -o <output file>]'
       sys.exit(2)
    print dataDir
    allContent, numberOfFeatures = readCSV(dataDir)
    nObjects = len(allContent) 

    subjectFeatures = computeAggregateFeatures(allContent,numberOfFeatures-1)

    feature_names = [ ["NumberOfPixels", "size_feature"], 
                      ["PhysicalSize", "size_feature"], 
                      ["NumberOfPixelsOnBorder", "size_feature"], 
                      ["FeretDiameter", "size_feature"], 
                      ["PrincipalMoments0", "shape_feature"], 
                      ["PrincipalMoments1", "shape_feature"], 
                      ["Elongation", "shape_feature"], 
                      ["Perimeter", "size_feature"], 
                      ["Roundness", "shape_feature"], 
                      ["EquivalentSphericalRadius", "size_feature"], 
                      ["EquivalentSphericalPerimeter", "size_feature"], 
                      ["EquivalentEllipsoidDiameter0", "size_feature"], 
                      ["EquivalentEllipsoidDiameter1","size_feature"], 
                      ["Flatness", "shape_feature"], 
                      ["meanR", "intensity_feature"], 
                      ["meanG", "intensity_feature"], 
                      ["meanB", "intensity_feature"], 
                      ["stdR", "intensity_feature"], 
                      ["stdG", "intensity_feature"], 
                      ["stdB", "intensity_feature"] ]

    dobj = {}
    dobj["patient_id"] = patientId 
    dobj["analysis_id"] = analysisId 
    dobj["visit_id"] = "visit-1"
    dobj["imaging_domain"] = "pathology"
    dobj["imaging_sequence"] = "H&E:tissue"
    nSubjectFeatures = len(subjectFeatures)    
    imaging_features = []
    i = 4
    for j in range(numberOfFeatures-5):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_Q25"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    i = i + 4
    for j in range(numberOfFeatures-5):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_Q75"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    i = i + 4
    for j in range(numberOfFeatures-5):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_median"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    dobj["imaging_features"] = imaging_features 

    fout = open(outFile,'w');
    fout.write(str(dobj))
    fout.write("\n");
    fout.close()

if __name__ == '__main__':
	main(sys.argv[1:])
