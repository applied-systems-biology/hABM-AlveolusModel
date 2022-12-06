add_library(visualisation SHARED
        PovFile.cpp
        PovRayObject.cpp
        Visualizer.cpp)
add_library(abm::visualisation ALIAS visualisation)
target_include_directories(visualisation PRIVATE ${PROJECT_SOURCE_DIR}/src)