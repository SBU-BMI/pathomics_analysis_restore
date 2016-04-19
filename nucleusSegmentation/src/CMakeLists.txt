cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(projName "imageAnalysis")
project(${projName})


set(CMAKE_C_FLAGS "-std=c99 -fPIC")
set(CMAKE_CXX_FLAGS " -fPIC")

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


################################################################################
# ITK, Should be with ItkVtkGlue
find_package (ITK)
if (ITK_FOUND)
   include(${ITK_USE_FILE})
endif(ITK_FOUND)


################################################################################
# put all libs together
set(Libraries
    ${ITK_LIBRARIES}
    ${Libraries}
)


################################################################################
# openCV
find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})
set(Libraries
    ${OpenCV_LIBS}
    ${Libraries}
)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package (GLIB)
if (HAVE_GLIB)
include_directories(${GLIB_INCLUDES})
set(Libraries
    ${GLIB_LIBRARIES}
    ${Libraries}
)
endif(HAVE_GLIB)


option(build_mainTileAndSegmentWSINucleiOutputFeature "build mainTileAndSegmentWSINucleiOutputFeature?" ON)
################################################################################
# find OpenSlide path
if(build_mainTileAndSegmentWSINucleiOutputFeature)
	find_path(OpenSlide_INSTALLATION_PATH include/openslide/openslide.h REQUIRED)
	if(NOT OpenSlide_INSTALLATION_PATH)
   		message(FATAL_ERROR "Please enter the installation path of OpenSlide, i.e., the one after --prefix= when configuring OpenSlide")
	endif(NOT OpenSlide_INSTALLATION_PATH)

	include_directories(${OpenSlide_INSTALLATION_PATH}/include/openslide/)

	set(OpenSlideLibs
  		${OpenSlide_INSTALLATION_PATH}/lib/libopenslide.so
  	)
endif(build_mainTileAndSegmentWSINucleiOutputFeature)


################################################################################
# put all libs together
set(Libraries
    ${OpenSlideLibs}
    ${Libraries}
)


################################################################################
# build local libraries
add_library(localLibraries 
  MultipleObjectFeatureAnalysisFilter.cxx
  SingleObjectFeatureAnalysisFilter.cxx
  BinaryMaskAnalysisFilter.cxx
  InputParameters.cxx
  Normalization.cpp
  PixelOperations.cpp
  Logger.cpp)

set(Libraries
    localLibraries
    ${Libraries}
)

option(build_mainSegmentSmallImage "build mainSegmentSmallImage?" ON)
if (build_mainSegmentSmallImage)
    set(cexx mainSegmentSmallImage)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainSegmentSmallImage)

# option(build_mainTileAndSegmentWSINuclei "build mainTileAndSegmentWSINuclei?" ON)
# if (build_mainTileAndSegmentWSINuclei)
#     set(cexx mainTileAndSegmentWSINuclei)
#     add_executable(${cexx} ${cexx}.cxx)
#     target_link_libraries(${cexx} ${Libraries})
# endif (build_mainTileAndSegmentWSINuclei)

option(build_mainTileAndSegmentWSINucleiOutputFeature "build mainTileAndSegmentWSINucleiOutputFeature?" ON)
if (build_mainTileAndSegmentWSINucleiOutputFeature)
    set(cexx mainTileAndSegmentWSINucleiOutputFeature)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainTileAndSegmentWSINucleiOutputFeature)

option(build_mainComputeShapeFeatures "build mainComputeShapeFeatures?" ON)
if (build_mainComputeShapeFeatures)
    set(cexx mainComputeShapeFeatures)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainComputeShapeFeatures)

option(build_mainComputeTileFeatures "build mainComputeTileFeatures?" ON)
if (build_mainComputeTileFeatures)
    set(cexx mainComputeTileFeatures)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainComputeTileFeatures)

option(build_mainSegmentFeatures "build mainSegmentFeatures?" ON)
if (build_mainSegmentFeatures)
    set(cexx mainSegmentFeatures)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainSegmentFeatures)

