cmake_minimum_required(VERSION 2.8)
 
project(TextonNoise LANGUAGES CXX)

find_package(ASTex)

add_executable(texton_cpu texton_cpu.h texton_cpu.cpp histogram.h main.cpp)
target_link_libraries(texton_cpu ASTex)

