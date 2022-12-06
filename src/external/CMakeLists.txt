add_library(xml_parser SHARED
        xmlParser/xmlParser.cpp)
add_library(external::xml_parser ALIAS xml_parser)
target_include_directories(xml_parser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/xmlParser)