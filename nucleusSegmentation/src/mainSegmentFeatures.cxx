#include <cstdio>
#include <iostream>
#include <string>
#include <omp.h>

#include <fstream>
#include <sstream>
#include <iostream>

#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

//itk
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkOpenCVImageBridge.h"

// openslide
#include "openslide.h"

// normalization
#include "Normalization.h"

// openCV
#include <opencv2/opencv.hpp>

#include "itkTypedefs.h"
#include "BinaryMaskAnalysisFilter.h"
#include "SFLSLocalChanVeseSegmentor2D.h"
#include "utilityIO.h"
#include "utilityTileAnalysis.h"
#include "utilityScalarImage.h"
#include "utilityWholeSlideProcessing.h"
#include "MultipleObjectFeatureAnalysisFilter.h"

#include "InputParameters.h"

const int _numOfFeatures = 25;
const std::string _featureNames[] = {
		"BoundingBoxTopLeftX",
        "BoundingBoxTopLeftY",
        "BoundingBoxBottomRightX",
        "BoundingBoxBottomRightY",
        "SizeInPixels",
        "PhysicalSize",
        "NumberOfPixelsOnBorder",
        "FeretDiameter",
        "PrincipalMoments0",
        "PrincipalMoments1",
        "Elongation",
        "Perimeter",
        "Roundness",
        "EquivalentSphericalRadius",
        "EquivalentSphericalPerimeter",
        "EquivalentEllipsoidDiameter0",
        "EquivalentEllipsoidDiameter1",
        "Flatness",
		"meanR",
    	"meanG",
    	"meanB",
    	"stdR",
    	"stdG",
    	"stdB",
		"Polygon"
};

typedef InputParameters AnalysisParameters;
int writeAnalysisParametersJSON(AnalysisParameters *analysisParams)  
{
	std::ostringstream oss;
	oss << analysisParams->outPrefix 
		<< "_mpp_" << analysisParams->mpp
		<< "_x" << analysisParams->topLeftX
		<< "_y" << analysisParams->topLeftY
		<< "-algmeta.json";
	std::ofstream outputMetadataFile(oss.str().c_str());

	outputMetadataFile	<< "{ "
						<< "\"inputType\" : \"" << analysisParams->inpType << "\", " 
						<< "\"otsuRatio\" : " << analysisParams->otsuRatio << ", "
						<< "\"curvatureWeight\" : " << analysisParams->curvatureWeight << ", "
						<< "\"sizeLowerThld\" : " << analysisParams->sizeLowerThld << ", "
						<< "\"sizeUpperThld\" : " << analysisParams->sizeUpperThld << ", "
						<< "\"msKernel\" : " << analysisParams->msKernel << ", "
						<< "\"levelsetNumberOfIteration\" : " << analysisParams->levelsetNumberOfIteration << ", "
						<< "\"topLeftX\" : " << analysisParams->topLeftX << ", "
						<< "\"topLeftY\" : " << analysisParams->topLeftY << ", "
						<< "\"sizeX\" :  " << analysisParams->sizeX << ", "
						<< "\"sizeY\" : " << analysisParams->sizeY << ", "
						<< "\"mpp\" : " << analysisParams->mpp << ", "
						<< "\"tileSize\" : " << analysisParams->tileSize << ", "
						<< "\"outPrefix\" : \"" << analysisParams->outPrefix << "\", "
						<< "\"outputLevel\" : \"" << analysisParams->outputLevel << "\", "
						<< "\"analysisId\" : \"" << analysisParams->analysisId << "\", "
						<< "\"analysisDesc\" : \"" << analysisParams->analysisDesc << "\""
						<< " }" << std::endl;
	outputMetadataFile.close();
	return 0;
}

int captureAnalysisParameters(AnalysisParameters *analysisParams, InputParameters *inpParams)
{
	analysisParams->inpType = inpParams->inpType;
	analysisParams->otsuRatio = inpParams->otsuRatio;
	analysisParams->curvatureWeight = inpParams->curvatureWeight;
	analysisParams->sizeLowerThld = inpParams->sizeLowerThld;
	analysisParams->sizeUpperThld = inpParams->sizeUpperThld;
	analysisParams->msKernel = inpParams->msKernel;
  	analysisParams->levelsetNumberOfIteration = inpParams->levelsetNumberOfIteration;
	analysisParams->topLeftX = inpParams->topLeftX;
	analysisParams->topLeftY = inpParams->topLeftY;
	analysisParams->sizeX = inpParams->sizeX;
	analysisParams->sizeY = inpParams->sizeY;
	analysisParams->mpp = inpParams->mpp;
	analysisParams->tileSize = inpParams->tileSize;	
	analysisParams->outputLevel = inpParams->outputLevel;
	analysisParams->inpFile = inpParams->inpFile;
	analysisParams->outPrefix = inpParams->outPrefix;
	analysisParams->analysisId = inpParams->analysisId;
	analysisParams->analysisDesc = inpParams->analysisDesc;
	return 0;
}

int writeFeatureCSV(std::string outPrefix, float mpp, int64_t topLeftX, int64_t topLeftY, std::vector< std::vector<FeatureValueType> > &features)
{
	std::ostringstream oss;
	oss << outPrefix 
		<< "_mpp_" << mpp 
		<< "_x" << topLeftX
		<< "_y" << topLeftY
		<< "-features.csv";
	std::ofstream outputFeatureFile(oss.str().c_str());
	int i;
	for (i=0;i<_numOfFeatures-1;i++) 
		outputFeatureFile<<_featureNames[i]<<",";
	outputFeatureFile<<_featureNames[i]<<std::endl;

	for (std::size_t iObject = 0; iObject < features.size(); ++iObject) {
		for (std::size_t iFeature = 0; iFeature < features[iObject].size(); ++iFeature)
			outputFeatureFile<<features[iObject][iFeature]<<",";
		outputFeatureFile<<std::endl<<std::flush;
	}
	outputFeatureFile.close();
}

int segmentWSI(InputParameters *inpParams)
{
	openslide_t *osr = openslide_open(inpParams->inpFile.c_str());
	if (osr==NULL) return 1;

	AnalysisParameters analysisParams;
	captureAnalysisParameters(&analysisParams,inpParams);

	analysisParams.mpp = ImagenomicAnalytics::WholeSlideProcessing::extractMPP<char>(osr);

	int32_t levelOfLargestSize = 0; // 0-th level is the largest
	int64_t largestW;
	int64_t largestH;
	{
		int64_t w[1];
		int64_t h[1];

		openslide_get_level_dimensions(osr, levelOfLargestSize, w, h);

		largestW = w[0];
		largestH = h[0];
	}

	std::vector<int64_t> tileTopleftX;
	std::vector<int64_t> tileTopleftY;
	std::vector<int64_t> tileSizeX;
	std::vector<int64_t> tileSizeY;

	ImagenomicAnalytics::WholeSlideProcessing::generateTileRegions<char>(largestW, largestH, inpParams->tileSize, 
													tileTopleftX, tileTopleftY, tileSizeX, tileSizeY);
	std::size_t numberOfTiles = tileSizeY.size();

#pragma omp parallel for
	for (std::size_t iTile = 0; iTile < numberOfTiles; ++iTile)
	{
		int64_t topLeftX = tileTopleftX[iTile];
		int64_t topLeftY = tileTopleftY[iTile];
		int64_t sizeX = tileSizeX[iTile];
		int64_t sizeY = tileSizeY[iTile];

		cv::Mat thisTile;
#pragma omp critical
		{
			thisTile = ImagenomicAnalytics::WholeSlideProcessing::extractTileFromWSI<char>(osr, levelOfLargestSize, topLeftX, topLeftY, sizeX, sizeY);
		}
		itkRGBImageType::Pointer thisTileItk =  itk::OpenCVImageBridge::CVMatToITKImage< itkRGBImageType >( thisTile );

		itkUShortImageType::Pointer outputLabelImageUShort = itkUShortImageType::New();
		itkUCharImageType::Pointer nucleusBinaryMask = ImagenomicAnalytics::TileAnalysis::processTile<char>(thisTile, outputLabelImageUShort, 
				inpParams->otsuRatio, inpParams->curvatureWeight, inpParams->sizeLowerThld, 
				inpParams->sizeUpperThld, inpParams->mpp, inpParams->msKernel, inpParams->levelsetNumberOfIteration);

#pragma omp critical
		{
			if (inpParams->outputLevel>=MASK_IMG) {
				std::ostringstream oss;
				oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-tile.jpg"; // Output tile
				ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(thisTileItk, oss.str().c_str(), 0);
			}
			if (inpParams->outputLevel>=MASK_ONLY) {
				std::ostringstream oss;
				oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-seg.png"; // Mask tile
				ImagenomicAnalytics::IO::writeImage<itkUCharImageType>(nucleusBinaryMask, oss.str().c_str(), 0);
			}
			if (inpParams->outputLevel>=MASK_IMG_OVERLAY) {
				std::ostringstream oss;
				oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-overlay.jpg"; 
				itkRGBImageType::Pointer overlay = ImagenomicAnalytics::ScalarImage::generateSegmentationOverlay<char>(thisTileItk, nucleusBinaryMask);
				ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(overlay, oss.str().c_str(), 0);
			}
		}

		// Compute features
  		ImagenomicAnalytics::MultipleObjectFeatureAnalysisFilter featureAnalyzer;
  		featureAnalyzer.setInputRGBImage(thisTileItk);
  		featureAnalyzer.setObjectBinaryMask(nucleusBinaryMask);
		featureAnalyzer.setTopLeft(topLeftX,topLeftY);
  		featureAnalyzer.update();

  		std::vector< std::vector<FeatureValueType> > features = featureAnalyzer.getFeatures();

#pragma omp critical
		{
			writeFeatureCSV(inpParams->outPrefix, inpParams->mpp, topLeftX, topLeftY, features);
		}
	}

#pragma omp barrier
	openslide_close(osr);
	writeAnalysisParametersJSON(&analysisParams);

	return 0;
}

int segmentImg(InputParameters *inpParams)
{
	const int ImageDimension = 2;

	AnalysisParameters analysisParams;
	captureAnalysisParameters(&analysisParams,inpParams);

	cv::Mat thisTile = imread(inpParams->inpFile.c_str());
	itkRGBImageType::Pointer thisTileItk =  itk::OpenCVImageBridge::CVMatToITKImage< itkRGBImageType >( thisTile );
	analysisParams.sizeX = (int64_t)thisTile.cols;
	analysisParams.sizeY = (int64_t)thisTile.rows;

	itkUShortImageType::Pointer outputLabelImage  = itkUShortImageType::New();
	itkUCharImageType::Pointer  nucleusBinaryMask = ImagenomicAnalytics::TileAnalysis::processTile<char>(thisTile, outputLabelImage, 
			inpParams->otsuRatio, inpParams->curvatureWeight, inpParams->sizeLowerThld, 
			inpParams->sizeUpperThld, inpParams->mpp, inpParams->msKernel,
			inpParams->levelsetNumberOfIteration);

	if (inpParams->outputLevel>=MASK_ONLY) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << inpParams->topLeftX << "_y" << inpParams->topLeftY << "-seg.png";
		ImagenomicAnalytics::IO::writeImage<itkUCharImageType>(nucleusBinaryMask, oss.str().c_str(), 0);
	}
	if (inpParams->outputLevel>=MASK_IMG) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << inpParams->topLeftX << "_y" << inpParams->topLeftY << "-tile.jpg"; 
		ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(thisTileItk, oss.str().c_str(), 0);
	}
	if (inpParams->outputLevel>=MASK_IMG_OVERLAY) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << inpParams->topLeftX << "_y" << inpParams->topLeftY << "-overlay.jpg"; 
		itkRGBImageType::Pointer overlay = ImagenomicAnalytics::ScalarImage::generateSegmentationOverlay<char>(thisTileItk, nucleusBinaryMask);
		ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(overlay, oss.str().c_str(), 0);
	}

	// Compute features
	ImagenomicAnalytics::MultipleObjectFeatureAnalysisFilter featureAnalyzer;
	featureAnalyzer.setInputRGBImage(thisTileItk);
	featureAnalyzer.setObjectBinaryMask(nucleusBinaryMask);
	featureAnalyzer.setTopLeft(0,0);
	featureAnalyzer.update();

	std::vector< std::vector<FeatureValueType> > features = featureAnalyzer.getFeatures();

	writeFeatureCSV(inpParams->outPrefix, inpParams->mpp, inpParams->topLeftX, inpParams->topLeftY, features);
	writeAnalysisParametersJSON(&analysisParams);

	return 0;
}

std::size_t readTileList(InputParameters *inpParams, std::vector<string>& inpFiles, std::vector<string>& outPrefixes, 
				std::vector<int64_t>& topLeftX, std::vector<int64_t>& topLeftY, 
				std::vector<int64_t>& sizeX, std::vector<int64_t>& sizeY)
{
	std::ifstream infile(inpParams->inpFile.c_str());

	std::string line;
	while (std::getline(infile,line)) {
		std::istringstream ss(line);
		std::string token;

		// Input image file
		std::getline(ss, token, ',');
		inpFiles.push_back(token);

		// Output file prefix	
		std::getline(ss, token, ',');
		outPrefixes.push_back(token);

		// top left X,Y
		std::getline(ss, token, ',');
		topLeftX.push_back(atoi(token.c_str()));
		std::getline(ss, token, ',');
		topLeftY.push_back(atoi(token.c_str()));

		// width and height
		std::getline(ss, token, ',');
		sizeX.push_back(atoi(token.c_str()));
		std::getline(ss, token, ',');
		sizeY.push_back(atoi(token.c_str()));
	}

	return (std::size_t) inpFiles.size();
}

int segmentTiles(InputParameters *inpParams)
{
	std::vector<string>  inpFileArray; 
	std::vector<string>  outPrefixArray; 
	std::vector<int64_t> topLeftXArray;
	std::vector<int64_t> topLeftYArray; 
	std::vector<int64_t> sizeXArray;
	std::vector<int64_t> sizeYArray;
	std::size_t          tileCount;

	tileCount = readTileList(inpParams,inpFileArray,outPrefixArray,
			topLeftXArray,topLeftYArray,sizeXArray,sizeYArray);

#pragma omp parallel for
	for (std::size_t iTile=0;iTile<tileCount;iTile++) {
		std::string fileName  = inpFileArray[iTile];
		std::string outPrefix = outPrefixArray[iTile];
		int64_t topLeftX      = topLeftXArray[iTile]; 
		int64_t topLeftY      = topLeftYArray[iTile];
		int64_t sizeX         = sizeXArray[iTile]; 
		int64_t sizeY         = sizeYArray[iTile];

		AnalysisParameters analysisParams;
		captureAnalysisParameters(&analysisParams,inpParams);
		analysisParams.topLeftX  = topLeftX;
		analysisParams.topLeftY  = topLeftY;
		analysisParams.sizeX     = sizeX;
		analysisParams.sizeY     = sizeY;
		analysisParams.outPrefix = outPrefix;

#pragma omp critical
		{
			std::cout 	<< "INPUT READING: " << fileName << " " << outPrefix << " "
						<< topLeftX << " " << topLeftY <<  " "
						<< sizeX << " " << sizeY << std::endl;
		}

		cv::Mat thisTile;
		int32_t levelOfLargestSize = 0; // 0-th level is the largest
		float   mpp;
		char    noErrors = 1;
#pragma omp critical
		{
			int64_t w[1],h[1];
			try {
				openslide_t *osr = openslide_open(fileName.c_str());
				mpp = ImagenomicAnalytics::WholeSlideProcessing::extractMPP<char>(osr);
				openslide_get_level_dimensions(osr, levelOfLargestSize, w, h);
				if ((topLeftX+sizeX)>w[0] || (topLeftY+sizeY)>h[0]) {
					openslide_close(osr);
					throw 1;
				}
				thisTile = ImagenomicAnalytics::WholeSlideProcessing::extractTileFromWSI<char>(osr, levelOfLargestSize, 
							topLeftX, topLeftY, sizeX, sizeY);
				openslide_close(osr);
			} catch (...) {
				std::cerr 	<< "ERROR: Requested tile ("  
							<< topLeftX <<"," << topLeftY << "," << topLeftX+sizeX << "," << topLeftY+sizeY
							<< ") is out of bounds (" 
							<< w[0] << "," << h[0] 
							<< ") in image: " << fileName << std::endl;
				noErrors = 0; 
			}
		}
		if (noErrors) {
			itkRGBImageType::Pointer thisTileItk = itk::OpenCVImageBridge::CVMatToITKImage< itkRGBImageType >( thisTile );
			analysisParams.mpp = mpp;

			itkUShortImageType::Pointer outputLabelImageUShort = itkUShortImageType::New();
			itkUCharImageType::Pointer nucleusBinaryMask = ImagenomicAnalytics::TileAnalysis::processTile<char>(thisTile, outputLabelImageUShort, 
					inpParams->otsuRatio, inpParams->curvatureWeight, inpParams->sizeLowerThld, 
					inpParams->sizeUpperThld, mpp, inpParams->msKernel, inpParams->levelsetNumberOfIteration);

#pragma omp critical
			{
				if (inpParams->outputLevel>=MASK_IMG) {
					std::ostringstream oss;
					oss << outPrefix << "_mpp_" << mpp << "_x" << topLeftX << "_y" << topLeftY << "-tile.jpg"; 
					ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(thisTileItk, oss.str().c_str(), 0);
				}
				if (inpParams->outputLevel>=MASK_ONLY) {
					std::ostringstream oss;
					oss << outPrefix << "_mpp_" << mpp << "_x" << topLeftX << "_y" << topLeftY << "-seg.png"; // Mask tile 
					ImagenomicAnalytics::IO::writeImage<itkUCharImageType>(nucleusBinaryMask, oss.str().c_str(), 0);
				}
				if (inpParams->outputLevel>=MASK_IMG_OVERLAY) {
					std::ostringstream oss;
					oss << outPrefix << "_mpp_" << mpp << "_x" << topLeftX << "_y" << topLeftY << "-overlay.jpg"; 
					itkRGBImageType::Pointer overlay = ImagenomicAnalytics::ScalarImage::generateSegmentationOverlay<char>(thisTileItk, nucleusBinaryMask);
					ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(overlay, oss.str().c_str(), 0);
				}
			}

			// Compute features
			ImagenomicAnalytics::MultipleObjectFeatureAnalysisFilter featureAnalyzer;
			featureAnalyzer.setInputRGBImage(thisTileItk);
			featureAnalyzer.setObjectBinaryMask(nucleusBinaryMask);
			featureAnalyzer.setTopLeft(topLeftX,topLeftY);
			featureAnalyzer.update();
			std::vector< std::vector<FeatureValueType> > features = featureAnalyzer.getFeatures();

#pragma omp critical
			{
				writeFeatureCSV(outPrefix, mpp, topLeftX, topLeftY, features);
				writeAnalysisParametersJSON(&analysisParams);
			}
		}
	}

#pragma omp barrier
	return 0;
}

int segmentSingleTile(InputParameters *inpParams)
{
	std::string fileName  = inpParams->inpFile; 
	std::string outPrefix = inpParams->outPrefix;
	int64_t topLeftX = inpParams->topLeftX;  
	int64_t topLeftY = inpParams->topLeftY; 
	int64_t sizeX = inpParams->sizeX; 
	int64_t sizeY = inpParams->sizeY;

	AnalysisParameters analysisParams;
	captureAnalysisParameters(&analysisParams,inpParams);

	openslide_t *osr = openslide_open(fileName.c_str());
	if (osr==NULL) return 1;

	analysisParams.mpp = ImagenomicAnalytics::WholeSlideProcessing::extractMPP<char>(osr);

	int32_t levelOfLargestSize = 0; // 0-th level is the largest
	int64_t w[1],h[1];
	openslide_get_level_dimensions(osr, levelOfLargestSize, w, h);
	if ((topLeftX+sizeX)>w[0] || (topLeftY+sizeY)>h[0]) {
		std::cerr 	<< "ERROR: Requested tile ("  
					<< topLeftX <<"," << topLeftY << "," << topLeftX+sizeX << "," << topLeftY+sizeY
					<< ") is out of bounds (" 
					<< w[0] << "," << h[0] 
					<< ") in image: " << fileName << std::endl;
		return 1;
	}

	cv::Mat thisTile;
	thisTile = ImagenomicAnalytics::WholeSlideProcessing::extractTileFromWSI<char>(osr, levelOfLargestSize, 
					topLeftX, topLeftY, sizeX, sizeY);
	openslide_close(osr);

	itkRGBImageType::Pointer thisTileItk =  itk::OpenCVImageBridge::CVMatToITKImage< itkRGBImageType >( thisTile );

	itkUShortImageType::Pointer outputLabelImageUShort = itkUShortImageType::New();
	itkUCharImageType::Pointer nucleusBinaryMask = ImagenomicAnalytics::TileAnalysis::processTile<char>(thisTile, outputLabelImageUShort, 
			inpParams->otsuRatio, inpParams->curvatureWeight, inpParams->sizeLowerThld, 
			inpParams->sizeUpperThld, inpParams->mpp, inpParams->msKernel, inpParams->levelsetNumberOfIteration);

	if (inpParams->outputLevel>=MASK_IMG) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-tile.jpg"; 
		ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(thisTileItk, oss.str().c_str(), 0);
	}
	if (inpParams->outputLevel>=MASK_ONLY) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-seg.png"; // Mask tile 
		ImagenomicAnalytics::IO::writeImage<itkUCharImageType>(nucleusBinaryMask, oss.str().c_str(), 0);
	}
	if (inpParams->outputLevel>=MASK_IMG_OVERLAY) {
		std::ostringstream oss;
		oss << inpParams->outPrefix << "_mpp_" << inpParams->mpp << "_x" << topLeftX << "_y" << topLeftY << "-overlay.jpg"; 
		itkRGBImageType::Pointer overlay = ImagenomicAnalytics::ScalarImage::generateSegmentationOverlay<char>(thisTileItk, nucleusBinaryMask);
		ImagenomicAnalytics::IO::writeImage<itkRGBImageType>(overlay, oss.str().c_str(), 0);
	}

	// Compute features
	ImagenomicAnalytics::MultipleObjectFeatureAnalysisFilter featureAnalyzer;
	featureAnalyzer.setInputRGBImage(thisTileItk);
	featureAnalyzer.setObjectBinaryMask(nucleusBinaryMask);
	featureAnalyzer.setTopLeft(topLeftX,topLeftY);
	featureAnalyzer.update();
	std::vector< std::vector<FeatureValueType> > features = featureAnalyzer.getFeatures();

	writeFeatureCSV(inpParams->outPrefix, inpParams->mpp, topLeftX, topLeftY, features);
	writeAnalysisParametersJSON(&analysisParams);

	return 0;
}

int main(int argc, char **argv)
{
	const int ImageDimension = 2;
	InputParameters inpParams;

	if (parseInputParameters(argc,argv,&inpParams)!=0) {
		printParseError(argv);
		exit(-1);
	}
	printInputParameters(&inpParams);

	if (inpParams.inpType==WSI) {
		segmentWSI(&inpParams);
	} else if (inpParams.inpType==IMG) {
		segmentImg(&inpParams);
	} else if (inpParams.inpType==TILES) {
		segmentTiles(&inpParams);
	} else if (inpParams.inpType==ONETILE) {
		segmentSingleTile(&inpParams);
	} else {
		std::cerr << "Unknown input type." << std::endl;
		printParseError(argv);
		exit(-1);
	}

	return 0;
}
