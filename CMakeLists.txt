cmake_minimum_required(VERSION 3.10)

project(MALVA VERSION 1.3.1)

file(GLOB SOURCES *.cpp)

add_executable(malva-geno ${SOURCES})

set_target_properties(malva-geno
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)

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
find_path(KMC_INCLUDE_DIR NAMES kmc_file.h PATHS ${CMAKE_SYSTEM_PREFIX_PATH} PATH_SUFFIXES "kmc")
if (NOT KMC_INCLUDE_DIR)
    message(FATAL_ERROR "kmc_file.h (part of libkmc-dev) required")
else ()
    message(STATUS "KMC includes found ${KMC_INCLUDE_DIR}")
endif ()

find_library(HTSLIB hts REQUIRED)
if (NOT HTSLIB)
    message(FATAL_ERROR "libhts (libhts-dev) required")
else ()
    message(STATUS "HTSLIB found ${HTSLIB}")
endif ()

target_include_directories(malva-geno
       PRIVATE ${PROJECT_SOURCE_DIR}
       PRIVATE ${KMC_INCLUDE_DIR}
)

target_link_libraries(malva-geno
        PUBLIC sdsl
        PUBLIC divsufsort
        PUBLIC divsufsort64
        PUBLIC kmc
        PUBLIC hts
        PUBLIC z
        PRIVATE OpenMP::OpenMP_CXX)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")