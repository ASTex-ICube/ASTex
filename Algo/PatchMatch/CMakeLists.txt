cmake_minimum_required(VERSION 2.8)
 
project(PatchMatch LANGUAGES CXX)

find_package(ASTex)

add_executable(test_patchmatch pm_minimal.cpp)
target_link_libraries(test_patchmatch ASTex)
