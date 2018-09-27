cmake_minimum_required(VERSION 3.0)
 
project(ASTexBench LANGUAGES CXX)

find_package(ASTex)

add_executable(bench_traversals "")
target_sources(bench_traversals PRIVATE bench_traversals.cpp)
target_link_libraries(bench_traversals PRIVATE ASTex)

# nice vs solution
set_target_properties( bench_traversals PROPERTIES FOLDER Bench)