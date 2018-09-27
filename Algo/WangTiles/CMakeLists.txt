cmake_minimum_required(VERSION 3.0)

project(wangTiles)

find_package(Qt5Widgets REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(ASTex REQUIRED)

add_executable(wang "") 
target_sources(wang PRIVATE wang.cpp min_path_cut.h imageviewer.h)
target_link_libraries(wang PRIVATE ASTex Qt5::Widgets)

# nice vs solution
set_target_properties(wang PROPERTIES FOLDER WangTiles)
