cmake_minimum_required(VERSION 3.5)
 
project(PatchExchange LANGUAGES CXX)

#find_package(ASTex)  ## done outside

add_executable(cyclostationary_noise "")

target_sources(cyclostationary_noise PRIVATE cyclostationary_noise.cpp)

target_link_libraries(cyclostationary_noise PRIVATE ASTex)

# nice vs solution
set_target_properties(cyclostationary_noise PROPERTIES FOLDER Algo)
