cmake_minimum_required(VERSION 2.8)

project(TerrainEnhancement LANGUAGES CXX)

find_package(ASTex)

add_executable(test_terrain_enhancement "")
target_sources(test_terrain_enhancement PRIVATE test_terrain_enhancement.cpp)
target_link_libraries(test_terrain_enhancement PRIVATE ASTex)

# nice vs solution
set_target_properties(test_terrain_enhancement PROPERTIES FOLDER Algo)

