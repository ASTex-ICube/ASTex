cmake_minimum_required(VERSION 2.8)
 
project(TextonNoiseAnalysis LANGUAGES CXX)

find_package(ASTex)

add_executable(test_analyse_texton)
target_sources(test_analyse_texton PRIVATE analyse_texton.cpp analyse_texton.h test_analyse_texton.cpp)
target_link_libraries(test_analyse_texton PRIVATE ASTex)

