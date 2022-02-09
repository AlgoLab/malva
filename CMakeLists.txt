cmake_minimum_required(VERSION 3.10)

project(MALVA VERSION 1.3.1)

add_executable(malva-geno zstdstream.cpp MurmurHash3.cpp main.cpp)

set_target_properties(malva-geno
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)

find_package(ZLIB)

find_library(SDSL libsdsl.a REQUIRED)
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

find_library(ZSTD libzstd.a REQUIRED)
if (NOT ZSTD)
    message(FATAL_ERROR "ZSTD (libzstd-dev) required")
else ()
    message(STATUS "ZSTD found ${ZSTD}")
endif ()

target_include_directories(malva-geno
       PRIVATE ${PROJECT_SOURCE_DIR}
       PRIVATE ${KMC_INCLUDE_DIR}
)

target_link_libraries(malva-geno
        PUBLIC ${SDSL}
        PUBLIC kmc
        PUBLIC hts
        PUBLIC z
        PUBLIC ${ZSTD}
        PUBLIC pthread)