cmake_minimum_required(VERSION 3.0)
 
project(PerioSmooth LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


add_executable(perio_smooth_color perio_smooth_color.cpp)
target_link_libraries(perio_smooth_color ASTex)

add_executable(perio_smooth perio_smooth.cpp)
target_link_libraries(perio_smooth ASTex)


# nice vs solution
set_target_properties(perio_smooth perio_smooth_color PROPERTIES FOLDER Perio_Smooth)


