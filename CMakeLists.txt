cmake_minimum_required(VERSION 3.15)
project(CEIM_EcoDeviceScore LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CEIM_ECODEVICESCORE_BUILD_EXAMPLES "Build example executable" ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if(CEIM_ECODEVICESCORE_BUILD_EXAMPLES)
    add_executable(ceim_ecodevice_score_example
        examples/ceim_ecodevice_score_example.cpp
    )
    target_compile_features(ceim_ecodevice_score_example PUBLIC cxx_std_17)
endif()
