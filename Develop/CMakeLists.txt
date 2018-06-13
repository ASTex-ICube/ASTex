cmake_minimum_required(VERSION 2.8)
 
project(Develop LANGUAGES CXX)

find_package(Qt5Widgets REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(ASTex)

add_executable(texton_cpu ContentExchange/content.h ContentExchange/patch.h ContentExchange/patchProcessor.h Stamping/stamper.h Stamping/stamp.h Stamping/sampler.h Stamping/sampler.cpp mipmap.h histogram.h rpn_utils.h rpn_utils.cpp texton_io.h texton_io.cpp gaussian_process.h imageviewer.h main.cpp)
#target_link_libraries(texton_cpu ASTex)
target_link_libraries(texton_cpu ASTex Qt5::Widgets)
