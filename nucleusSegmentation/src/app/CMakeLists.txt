option(build_mainSegmentSmallImage "build mainSegmentSmallImage?" ON)
if (build_mainSegmentSmallImage)
    set(cexx mainSegmentSmallImage)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainSegmentSmallImage)

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

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

option(build_mainSegmentFeatures "build mainSegmentFeatures?" ON)
if (build_mainSegmentFeatures)

	ExternalProject_Add(nscale-proj
		DOWNLOAD_DIR nscale
		GIT_REPOSITORY https://github.com/SBU-BMI/nscale.git 
		SOURCE_DIR nscale
		CMAKE_ARGS -DOpenCV_DIR=${OpenCV_DIR}
		INSTALL_COMMAND ""
	)

	set(nscalesrcdir   ${CMAKE_BINARY_DIR}/app/nscale/src)
	set(nscalebuilddir ${CMAKE_BINARY_DIR}/app/nscale-proj-prefix/src/nscale-proj-build)

    #includes
    include_directories(${nscalesrcdir}/segment)
    include_directories(${nscalesrcdir}/segment/cuda)
    include_directories(${nscalesrcdir}/common)
    include_directories(${nscalesrcdir}/features)
    include_directories(${nscalesrcdir}/normalization)
	include_directories(../lib/ConvertFunctions)

    link_directories(${nscalebuilddir}/lib)

    set(cexx mainSegmentFeatures)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} segment featuresAll ${Libraries})
	add_dependencies(mainSegmentFeatures nscale-proj)
endif (build_mainSegmentFeatures)

# if (build_mainSegmentFeatures)
# 	add_subdirectory(ConvertFunctions)
# endif (build_mainSegmentFeatures)


