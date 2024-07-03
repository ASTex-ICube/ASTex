cmake_minimum_required(VERSION 3.5)

project(wangTiles)

find_package(ASTex REQUIRED)
find_package(Qt5 OPTIONAL_COMPONENTS Widgets QUIET)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


if (Qt5Widgets_FOUND)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  add_executable(wang_qt "")
  target_sources(wang_qt PUBLIC wang_qt.cpp wang_tiles.h min_path_cut.h imageviewer.h)
  target_link_libraries(wang_qt Qt5::Widgets ASTex)
  set_target_properties(wang_qt PROPERTIES FOLDER Algo)
endif()

add_executable(wang "")
target_sources(wang PRIVATE wang.cpp wang_tiles.h min_path_cut.h)
target_link_libraries(wang ASTex)


# nice vs solution
set_target_properties(wang PROPERTIES FOLDER Algo)
