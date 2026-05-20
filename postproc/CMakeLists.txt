# Minimal CMake file to compile cross_plane.cpp file

cmake_minimum_required(VERSION 3.10)
project(cross_plane LANGUAGES CXX)
add_executable(cross_plane cross_plane.cpp)
set_target_properties(cross_plane PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
install(TARGETS cross_plane
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
