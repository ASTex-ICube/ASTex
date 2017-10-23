cmake_minimum_required(VERSION 3.0)
 
project(Quilting LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


add_executable(test_quilting test_quilting.cpp quilting.cpp quilting.h error_buffer.cpp error_buffer.h)
target_link_libraries(test_quilting ASTex)

# nice vs solution
set_target_properties(test_quilting PROPERTIES FOLDER Quilting)

