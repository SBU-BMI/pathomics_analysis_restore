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

option(build_mainSegmentFeatures "build mainSegmentFeatures?" ON)
if (build_mainSegmentFeatures)
    set(cexx mainSegmentFeatures)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} ${Libraries} )
endif (build_mainSegmentFeatures)

