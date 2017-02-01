set(SOURCE_FILES
        ConvertBinaryToLabel.cpp
        ConvertLabelToCvMat.cpp
        ConvertRGBToLabel.cpp
        ConvertSlicerToLabel.cpp
        IsImageBinary.cpp
        ReadUnknownImageType.cpp
        )

include_directories("../../include")
	
add_library(ConvertFunctions ${SOURCE_FILES})


