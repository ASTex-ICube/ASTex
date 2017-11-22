cmake_minimum_required(VERSION 3.0)

project(wangTiles)

find_package(ASTex REQUIRED)
find_package(Qt5Widgets)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )



if (Qt5Widgets_FOUND)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  add_executable(wang_qt wang_qt.cpp wang_tiles.h min_path_cut.h imageviewer.h)
  target_link_libraries(wang_qt Qt5::Widgets ASTex)
endif()

add_executable(wang wang.cpp wang_tiles.h min_path_cut.h imageviewer.h)
target_link_libraries(wang Qt5::Widgets ASTex)


