cmake_minimum_required(VERSION 3.0)

project(ProceduralNoiseFiltering)

find_package(ASTex REQUIRED)
find_package(Qt5 COMPONENTS Widgets)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


if (Qt5Widgets_FOUND)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  add_executable(pnf_qt "")
  target_sources(pnf_qt PUBLIC pnf_qt.cpp imageviewer.h texture_noise.h pnf.h)
  target_link_libraries(pnf_qt Qt5::Widgets ASTex)
endif()

add_executable(pnf "")
target_sources(pnf PRIVATE pnf.cpp pnf.h)
target_link_libraries(pnf ASTex)


# nice vs solution
set_target_properties(pnf PROPERTIES FOLDER ProceduralNoiseFiltering)
