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
       f = open(j,'r')
       readContent = f.read().split('\n')
       numFeatures = len(readContent[0].split(','))
       print numFeatures
       featureList = readContent[0].split(',')
       if featureList[numFeatures-1]=='' or featureList[numFeatures-1]==' ':
          numFeatures = numFeatures-1
       readContent.pop(0)
       allContent += readContent 
       mydata = allContent.pop()
       f.close()   
    return (allContent,numFeatures)

def writeHeaderCSV(feature_names,fout): 
    fout.write("ParticipantBarcode")
    fout.write(",")
    fout.write("ImagingModality")
    fout.write(",")
    featureLen = len(feature_names)
    for j in range(featureLen):
        fout.write(str(feature_names[j][0]) + "_Q25")
        fout.write(",")
    for j in range(featureLen):
        fout.write(str(feature_names[j][0]) + "_median")
        fout.write(",")
    for j in range(featureLen):
        fout.write(str(feature_names[j][0]) + "_Q75")
        fout.write(",")
    fout.write("SampleTypeLetterCode,")
    fout.write("Study,")
    fout.write("AnalysisId\n")

def writeDataCSV(feature_names,patientId,subjectFeatures,numberOfFeatures,skipBegin,skipEnd,studyId,letterCode,analysisId,fout):
    writeHeaderCSV(feature_names,fout)
    fout.write(str(patientId))
    fout.write(",")
    fout.write("pathology")
    fout.write(",")
    nSubjectFeatures = len(subjectFeatures)    
    featureLen = numberOfFeatures-(skipEnd+skipBegin);
    i = skipBegin
    for j in range(featureLen):
        fout.write(str(subjectFeatures[i]))
        fout.write(",")
        i = i + 1
    i = i + skipBegin
    for j in range(featureLen):
        fout.write(str(subjectFeatures[i]))
        fout.write(",")
        i = i + 1
    i = i + skipBegin
    for j in range(featureLen):
        fout.write(str(subjectFeatures[i]))
        fout.write(",")
        i = i + 1
    fout.write(str(letterCode))
    fout.write(",")
    fout.write(str(studyId))
    fout.write(",")
    fout.write(str(analysisId))
    fout.write("\n")

def writeDataJSON(feature_names,patientId,subjectFeatures,numberOfFeatures,skipBegin,skipEnd,studyId,letterCode,analysisId,fout):
    dobj = {}
    dobj["patient_id"] = patientId 
    dobj["analysis_id"] = analysisId 
    dobj["visit_id"] = "visit-1"
    dobj["imaging_domain"] = "pathology"
    dobj["imaging_sequence"] = "H&E:tissue"

    nSubjectFeatures = len(subjectFeatures)    
    imaging_features = []
    featureLen = numberOfFeatures-(skipEnd+skipBegin);
    print featureLen
    print len(feature_names)
    i = skipBegin
    for j in range(featureLen):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_Q25"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    i = i + skipBegin
    for j in range(featureLen):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_median"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    i = i + skipBegin
    for j in range(featureLen):
        dobj2 = {}
        dobj2["feature_name"] = feature_names[j][0] + "_Q75"
        dobj2["value"] = subjectFeatures[i] 
        dobj2["feature_type"] = feature_names[j][1]
        imaging_features.append(dobj2)
        i = i + 1
    dobj["imaging_features"] = imaging_features 
    dobj["cancer_type"] = studyId
    dobj["tumor_type"]  = letterCode

    fout.write(str(dobj))
    fout.write("\n");

def computeAggregateFeatures(allContent,numberOfFeatures):
    nObjects = len(allContent)
    allFeatureMatrix = np.empty([nObjects, numberOfFeatures], dtype=float)

    for iObject in range(nObjects):
        featuresOfThisObject = allContent[iObject]
        featuresOfThisObject = featuresOfThisObject.split(',')
        featuresOfThisObject = featuresOfThisObject[:numberOfFeatures]
        featuresOfThisObject = [float(x) for x in featuresOfThisObject]

        allFeatureMatrix[iObject, :] = np.asarray(featuresOfThisObject)

    [featureQ25, featureMedian, featureQ75] = np.percentile(allFeatureMatrix, [25, 50, 75], axis=0)
    subjectFeatures = np.concatenate((featureQ25, featureMedian, featureQ75))
    return subjectFeatures

def printArgHelp():
    print 'mainAggregateFeatures.py [-h] [-i <inputfile pattern> -p <patient id> -a <analysis id> -o <output file> -t <out type: csv|json> -s <cancer type> -l <tumor type>]'
    print '\n'

def main(argv):
    dataDir = ''
    patientId = ''
    analysisId = ''
    outFile = ''
    cancerType = ''
    tumorType = ''
    outType = 'json'
    try:
		opts, args = getopt.getopt(argv,"hi:p:a:o:t:s:l:")
    except getopt.GetoptError:
       printArgHelp()
       sys.exit(2)
    if len(argv) == 0:
       printArgHelp()
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          printArgHelp()
          sys.exit()
       elif opt == "-i": 
          dataDir = arg
       elif opt == '-p':
          patientId = arg
       elif opt == '-a':
          analysisId = arg
       elif opt == '-o':
          outFile = arg
       elif opt == '-s':
          cancerType = arg
       elif opt == '-l':
          tumorType = arg
       elif opt == '-t':
          outType = arg

    if dataDir=='' or patientId=='' or analysisId=='' or cancerType=='' or outFile=='':
       printArgHelp()
       sys.exit(2)

    # Specific to mainSegmentFeatures output
    skipBegin = 4    # bounding box columns
    skipEnd   = 1    # polygon column

    feature_names = [ ["SizeInPixels", "size_feature"], 
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

    # Process the files
    allContent, numberOfFeatures = readCSV(dataDir)
    nObjects = len(allContent) 
    subjectFeatures = computeAggregateFeatures(allContent,numberOfFeatures-skipEnd)

    fout = open(outFile,'w');
    print outType
    if outType=='json':
       writeDataJSON(feature_names,patientId,subjectFeatures,numberOfFeatures,skipBegin,skipEnd,cancerType,tumorType,analysisId,fout)
    else:
       writeDataCSV(feature_names,patientId,subjectFeatures,numberOfFeatures,skipBegin,skipEnd,cancerType,tumorType,analysisId,fout)
    fout.close()

if __name__ == '__main__':
	main(sys.argv[1:])
