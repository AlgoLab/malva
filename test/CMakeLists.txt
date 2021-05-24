find_package(Catch2)

find_package(OpenMP)
find_package(ZLIB)

if (Catch2_FOUND)
    message(STATUS "Found Catch2, testing is enabled")
    enable_testing()
    file(GLOB TEST_SRC *.cpp)
    set(MALVA_SRC ${SOURCES})
    list(REMOVE_ITEM MALVA_SRC ${CMAKE_SOURCE_DIR}/main.cpp)
    foreach (TEST ${TEST_SRC})
        get_filename_component(TEST_NAME ${TEST} NAME_WE)
        add_executable(${TEST_NAME} ${TEST} ${MALVA_SRC})
        target_link_libraries(${TEST_NAME}
                PUBLIC Catch2::Catch2
                PUBLIC sdsl
                PUBLIC divsufsort
                PUBLIC divsufsort64
                PUBLIC kmc
                PUBLIC hts
                PUBLIC z
                PRIVATE OpenMP::OpenMP_CXX)
    endforeach ()
else ()
    message(WARNING "Catch2 not found, testing is disabled")
endif ()
