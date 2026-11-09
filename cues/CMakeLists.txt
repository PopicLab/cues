cmake_minimum_required(VERSION 3.16)
project(cues)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

find_package(ZLIB REQUIRED)
find_path(MINIMAP2_INCLUDE_DIR minimap.h PATHS $ENV{CONDA_PREFIX}/include)
find_library(MINIMAP2_LIBRARY minimap2 PATHS $ENV{CONDA_PREFIX}/lib)
add_library(minimap2 STATIC IMPORTED)
set_target_properties(minimap2 PROPERTIES
        IMPORTED_LOCATION ${MINIMAP2_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${MINIMAP2_INCLUDE_DIR}
)

find_path(HTSLIB_INCLUDE_DIR htslib/faidx.h PATHS $ENV{CONDA_PREFIX}/include)
find_library(HTSLIB_LIBRARY hts PATHS $ENV{CONDA_PREFIX}/lib)

######################
####### OpenMP #######
######################
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# C++ library
include_directories($ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(SOURCES src/cues.cc src/signals.cc)
add_library(${PROJECT_NAME} STATIC ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(${PROJECT_NAME} PUBLIC ${Python_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE minimap2 ZLIB::ZLIB)
target_include_directories(${PROJECT_NAME} PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${HTSLIB_LIBRARY})

add_executable (${PROJECT_NAME}_example src/main.cc)
target_link_libraries (${PROJECT_NAME}_example PUBLIC ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}_example PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Python wrapper
pybind11_add_module(py${PROJECT_NAME} src/bindings.cc)
target_link_libraries(py${PROJECT_NAME} PRIVATE ${PROJECT_NAME})
install(TARGETS py${PROJECT_NAME} LIBRARY DESTINATION .)


if(CUES_BUILD_TESTS)
    message(STATUS "Building tests...")
    enable_testing()
    add_subdirectory(tests)
endif()
