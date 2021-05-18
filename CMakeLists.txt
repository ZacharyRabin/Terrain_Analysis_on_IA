cmake_minimum_required(VERSION 3.1)

project(LibTri-Terrain-Analysis)
set(LibTri-Terrain-Analysis_VERSION_MAJOR 1)
set(LibTri-Terrain-Analysis_VERSION_MINOR 0)

set(LANGUAGE C++)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
    message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O3") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s -lrt") ## Strip binary
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-O3 -march=native") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s") ## Strip binary
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "/O2 /EHsc")
endif()

include_directories(sources)
add_subdirectory(sources)

### MAIN OF THE LIBRARY
set(EXECUTABLE_OUTPUT_PATH bin)
add_executable(ia_terrain_analysis sources/main.cpp sources/topological_main.cpp)
target_link_libraries(ia_terrain_analysis _terrain_analysis)

install (TARGETS ia_terrain_analysis RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})



