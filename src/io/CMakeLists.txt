include_directories(${Boost_INCLUDE_DIRS})
add_library(io SHARED
        InputConfiguration.cpp
        output_handler.cpp
        XMLFile.cpp)
add_library(abm::io ALIAS io)
target_include_directories(io PRIVATE ${PROJECT_SOURCE_DIR}/src)