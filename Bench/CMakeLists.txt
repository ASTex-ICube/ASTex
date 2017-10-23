cmake_minimum_required(VERSION 3.0)
 
project(ASTexBench LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\))


add_executable(bench_traversals bench_traversals.cpp)
target_link_libraries(bench_traversals ASTex)

# nice vs solution
set_target_properties( bench_traversals PROPERTIES FOLDER Bench)