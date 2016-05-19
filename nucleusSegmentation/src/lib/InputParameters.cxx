#include <cstdio>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "InputParameters.h"

void printParseError(char *argv[]) 
{
	std::cerr 	<< "Usage: " << argv[0] << " [parameters]" << std::endl;
	std::cerr	<< "   -t [wsi|tiles|onetile|img] " << std::endl
				<< "   -i <input file>" << std::endl
				<< "   -p <output file/prefix>" << std::endl
				<< "   -m <mpp>" << std::endl
				<< "   -r <otsuRatio> " << std::endl
				<< "   -w <curvatureWeight>" << std::endl
				<< "   -l <sizeLowerThld>" << std::endl
				<< "   -u <sizeUpperThld>" << std::endl
				<< "   -k <msKernel>" << std::endl
				<< "   -n <levelsetNumberOfIterations>" << std::endl
				<< "   -s <topLeftX,topLeftY>" << std::endl
				<< "   -b <sizeX,sizeY>" << std::endl
				<< "   -a <analysisId: string>" << std::endl
				<< "   -e <analysis desc: string>" << std::endl
				<< "   -d <tileSizeX,tileSizeY>" << std::endl
				<< "   -z <zipFile - compression works with onetile only.>" << std::endl
				<< "   -v <output level: mask|mask:img|mask:img:overlay>" << std::endl;
}

void printInputParameters(InputParameters *inpParams)
{
	std::cout << "Input type: " << inpParams->inpType << std::endl; 
	std::cout << "otsuRatio : " << inpParams->otsuRatio << std::endl;
	std::cout << "curvatureWeight: " << inpParams->curvatureWeight << std::endl;
	std::cout << "sizeLowerThld: " << inpParams->sizeLowerThld << std::endl;
	std::cout << "sizeUpperThld: " << inpParams->sizeUpperThld << std::endl;
	std::cout << "msKernel: " << inpParams->msKernel << std::endl;
  	std::cout << "levelsetNumberOfIteration: " << inpParams->levelsetNumberOfIteration << std::endl;
	std::cout << "topLeftX: " << inpParams->topLeftX << " topLeftY: " <<  inpParams->topLeftY << std::endl;
	std::cout << "sizeX:  " << inpParams->sizeX << " sizeY: " << inpParams->sizeY << std::endl;
	std::cout << "mpp: " << inpParams->mpp << std::endl;
	std::cout << "tileSizeX: " << inpParams->tileSizeX << " tileSizeY: " << inpParams->tileSizeY << std::endl;	
	std::cout << "inpFile: " << inpParams->inpFile << std::endl;
	std::cout << "outPrefix: " << inpParams->outPrefix << std::endl;
	std::cout << "outputLevel: " << inpParams->outputLevel << std::endl;
	if (inpParams->isZipped) std::cout << "zipFile: " << inpParams->zipFile;
}

int parseInputParameters(int argc, char **argv, InputParameters *inpParams) 
{
	int  c;

	if (argc<4) return 1;

	inpParams->inpType = WSI;
	inpParams->otsuRatio = 1.0;
	inpParams->curvatureWeight = 0.8;
	inpParams->sizeLowerThld = 3.0;
	inpParams->sizeUpperThld = 200.0;
	inpParams->msKernel = 20.0;
	inpParams->levelsetNumberOfIteration = 100;
	inpParams->topLeftX = 0;
	inpParams->topLeftY = 0;
	inpParams->sizeX = DEFAULT_SMALL_TILE;
	inpParams->sizeY = DEFAULT_SMALL_TILE;
	inpParams->mpp   = 0.25; // 40x objective
	inpParams->tileSizeX = 0;
	inpParams->tileSizeY = 0;
	inpParams->outputLevel = MASK_ONLY;
	inpParams->isZipped = 0; // no compressed zip output

	opterr = 0;
	while ((c = getopt (argc, argv, "ht:i:p:m:r:w:l:u:k:n:s:b:d:v:a:e:z:")) != -1) {
		switch (c)
		{
			case 'h': 
				return 1;	
			case 't': {
				if (!strcmp(optarg,"wsi")) {
					inpParams->inpType   = WSI;
				} else if (!strcmp(optarg,"tiles")) {
					inpParams->inpType   = TILES;
				} else if (!strcmp(optarg,"onetile")) {
					inpParams->inpType   = ONETILE;
				} else if (!strcmp(optarg,"img")) {
					inpParams->inpType = IMG;
				} else {
					fprintf(stderr, "Undefined input type.\n");
					return 1;
				}
				break;
			}
			case 'i':
				inpParams->inpFile = optarg;
				break;
			case 'p':
				inpParams->outPrefix = optarg;
				break;
			case 'z':
				inpParams->isZipped = 1;
				inpParams->zipFile = optarg;
				break;
			case 'a':
				inpParams->analysisId = optarg;
				break;
			case 'e':
				inpParams->analysisDesc = optarg;
				break;
			case 'm':
				inpParams->mpp = atof(optarg);
				break;
			case 'r':
				inpParams->otsuRatio = atof(optarg);
				break;
			case 'w':
				inpParams->curvatureWeight = atof(optarg);
				break;
			case 'l':
				inpParams->sizeLowerThld = atof(optarg);
				break;
			case 'u':
				inpParams->sizeUpperThld = atof(optarg);
				break;
			case 'k':
				inpParams->msKernel = atof(optarg);
				break;
			case 'n':
				inpParams->levelsetNumberOfIteration = (int64_t) atoi(optarg);
				break;
			case 's': {
				std::istringstream ss(optarg);
				std::string token;
				if (std::getline(ss, token, ',')) {
					inpParams->topLeftX = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -s is missing <leftX,leftY> value.\n");
					return 1;
				}
				if (std::getline(ss, token, ',')) {
					inpParams->topLeftY = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -s is missing <leftX,leftY> value.\n");
					return 1;
				}
				break;
			}
			case 'b': {
				std::istringstream ss(optarg);
				std::string token;
				if (std::getline(ss, token, ',')) {
					inpParams->sizeX = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -b is missing <sizeX,sizeY> value.\n");
					return 1;
				}
				if (std::getline(ss, token, ',')) {
					inpParams->sizeY = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -b is missing <sizeX,sizeY> value.\n");
					return 1;
				}
				break;
			}
			case 'd': {
				std::istringstream ss(optarg);
				std::string token;
				if (std::getline(ss, token, ',')) {
					inpParams->tileSizeX = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -d is missing <tileSizeX,tileSizeY> value.\n");
					return 1;
				}
				if (std::getline(ss, token, ',')) {
					inpParams->tileSizeY = atoi(token.c_str());
				} else {
					fprintf(stderr,"ERROR: Option -d is missing <tileSizeX,tileSizeY> value.\n");
					return 1;
				}
				break;
			}
			case 'v':
				if (!strcmp(optarg,"mask")) {
					inpParams->outputLevel = MASK_ONLY;
				} else if (!strcmp(optarg,"mask:img")) {
					inpParams->outputLevel = MASK_IMG;
				} else if (!strcmp(optarg,"img:mask")) {
					inpParams->outputLevel = MASK_IMG;
				} else if (!strcmp(optarg,"mask:img:overlay")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else if (!strcmp(optarg,"mask:overlay:img")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else if (!strcmp(optarg,"img:mask:overlay")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else if (!strcmp(optarg,"img:overlay:mask")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else if (!strcmp(optarg,"overlay:mask:img")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else if (!strcmp(optarg,"overlay:img:mask")) {
					inpParams->outputLevel = MASK_IMG_OVERLAY;
				} else {
					fprintf(stderr, "Undefined output level value.\n");
					return 1;
				}
				break;
			default:
				return 1;	
		}
	}

	if (inpParams->tileSizeX==0 || inpParams->tileSizeY==0) {
		if (inpParams->inpType==WSI) {
			inpParams->tileSizeX = DEFAULT_WSI_TILE;
			inpParams->tileSizeY = DEFAULT_WSI_TILE;
		} else {
			inpParams->tileSizeX = DEFAULT_SMALL_TILE;
			inpParams->tileSizeY = DEFAULT_SMALL_TILE;
		}
	}

	if (inpParams->topLeftX<0 || inpParams->topLeftY<0 
		|| inpParams->sizeX<0 || inpParams->sizeY<0) {
		fprintf(stderr,"Error in input parameter values. Check -s and -b values; one or more values are negative.\n");
		return 1;
	}	

	if (opterr) return 1;

	return 0;
}
