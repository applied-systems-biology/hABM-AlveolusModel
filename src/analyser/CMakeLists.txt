include_directories(${Boost_INCLUDE_DIRS})
add_library(analyser SHARED
        histogram_measurement.cpp
        pair_measurment.cpp
        Analyser.cpp
        InSituMeasurements.cpp)
add_library(abm::analyser ALIAS analyser)
target_include_directories(analyser PRIVATE ${PROJECT_SOURCE_DIR}/src)