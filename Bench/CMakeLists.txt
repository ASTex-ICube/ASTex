cmake_minimum_required(VERSION 3.5)
 
project(ASTexBench LANGUAGES CXX)

find_package(ASTex)

add_executable(bench_traversals "")
target_sources(bench_traversals PRIVATE bench_traversals.cpp)
target_link_libraries(bench_traversals PRIVATE ASTex)

add_executable(bench_eigen_pixels bench_eigen_pixels.cpp)
target_link_libraries(bench_eigen_pixels ASTex)


#add_executable(bench_parallel_for bench_parallel_for.cpp)
#target_link_libraries(bench_parallel_for ASTex)

# nice vs solution
#set_target_properties( bench_traversals bench_parallel_for PROPERTIES FOLDER Bench)
set_target_properties( bench_traversals bench_eigen_pixels PROPERTIES FOLDER Bench)

