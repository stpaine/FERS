cmake_minimum_required(VERSION 3.0)
project(your_project_name)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

find_package(XercesC REQUIRED)
include_directories(${XercesC_INCLUDE_DIRS})
link_directories(${XercesC_LIBRARY_DIRS})

set(SRCS validator.cpp xml_validator_output.cpp kml_visualiser.cpp)

add_executable(validator validator.cpp)
target_link_libraries(validator ${XercesC_LIBRARIES})

add_executable(xml_validator_output xml_validator_output.cpp)
target_link_libraries(xml_validator_output ${XercesC_LIBRARIES})

add_executable(kml_visualiser kml_visualiser.cpp)
target_link_libraries(kml_visualiser ${XercesC_LIBRARIES})

install(TARGETS validator xml_validator_output kml_visualiser DESTINATION bin)
