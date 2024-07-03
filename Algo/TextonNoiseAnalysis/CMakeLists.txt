cmake_minimum_required(VERSION 3.5)
 
project(TextonNoiseAnalysis LANGUAGES CXX)

find_package(ASTex)

add_executable(test_analyse_texton "")
target_sources(test_analyse_texton PRIVATE analyse_texton.cpp analyse_texton.h test_analyse_texton.cpp)
target_link_libraries(test_analyse_texton PRIVATE ASTex)

# nice vs solution
set_target_properties(test_analyse_texton PROPERTIES FOLDER Algo)

