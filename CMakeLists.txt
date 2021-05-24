cmake_minimum_required(VERSION 3.10)

project(MALVA VERSION 1.3.1)

file(GLOB SOURCES *.cpp)

add_executable(malva-geno ${SOURCES})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

find_package(OpenMP)
find_package(ZLIB)

find_library(SDSL sdsl REQUIRED)
if (NOT SDSL)
    message(FATAL_ERROR "sdsl-lite (libsdsl-dev) required")
else ()
    message(STATUS "SDSL found ${SDSL}")
endif ()

find_library(KMC kmc REQUIRED)
if (NOT KMC)
    message(FATAL_ERROR "kmc (libkmc-dev) required")
else ()
    message(STATUS "KMC found ${KMC}")
endif ()

find_library(HTSLIB hts REQUIRED)
if (NOT HTSLIB)
    message(FATAL_ERROR "libhts (libhts-dev) required")
else ()
    message(STATUS "HTSLIB found ${HTSLIB}")
endif ()

include_directories("${PROJECT_SOURCE_DIR}")

target_link_libraries(malva-geno
        PUBLIC sdsl
        PUBLIC divsufsort
        PUBLIC divsufsort64
        PUBLIC kmc
        PUBLIC hts
        PUBLIC z
        PRIVATE OpenMP::OpenMP_CXX)
