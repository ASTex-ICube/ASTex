cmake_minimum_required(VERSION 2.8)
 
project(PatchMatch LANGUAGES CXX)

find_package(ASTex)

add_executable(test_patchmatch "")
target_sources(test_patchmatch PRIVATE pm_minimal.cpp)
target_link_libraries(test_patchmatch PRIVATE ASTex)

# nice vs solution
set_target_properties(test_patchmatch PROPERTIES FOLDER tests)