cmake_minimum_required(VERSION 3.13)
project(bhut)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)


include_directories("${CMAKE_CURRENT_SOURCE_DIR}/includes")

add_executable(bhut main.cpp)
target_compile_options(bhut PRIVATE -static-libstdc++ -static -O0 -g -pedantic -Wall)

add_executable(transform transform.cpp)

