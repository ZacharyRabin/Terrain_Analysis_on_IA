cmake_minimum_required(VERSION 3.1)

include_directories(sources)
add_subdirectory(sources)

target_include_directories(_IA2D PUBLIC sources)
