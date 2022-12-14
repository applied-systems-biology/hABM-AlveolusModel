add_subdirectory(analyser)
add_subdirectory(basic)
add_subdirectory(io)
add_subdirectory(simulation)
add_subdirectory(utils)
add_subdirectory(external)
add_subdirectory(visualisation)

target_link_libraries(analyser PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(basic PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(io PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(simulation PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(utils PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(visualisation PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)

add_executable(hABM main.cpp)
target_link_libraries(hABM PUBLIC
        project_options
        project_warnings
        abm::simulation
        abm::analyser
        abm::basic
        abm::io
        abm::utils
        abm::visualisation
        external::xml_parser
        Boost::filesystem
        OpenMP::OpenMP_CXX)

target_include_directories(hABM PRIVATE ${PROJECT_SOURCE_DIR}/src)


