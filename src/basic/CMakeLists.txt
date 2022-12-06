add_library(basic SHARED
        ColorRGB.cpp
        Coordinate3D.cpp
        InversionSampler.cpp
        Randomizer.cpp
        Sampler.cpp
        SphericCoordinate3D.cpp
        )
add_library(abm::basic ALIAS basic)
target_include_directories(basic PRIVATE ${PROJECT_SOURCE_DIR}/src)