cmake_minimum_required(VERSION 3.0)
 
project(Graphcut)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


find_package(ASTex)

add_executable(test_graphcut test_graphcut.cpp graphcut.h graphcut.cpp maxflow/block.h maxflow/graph.h maxflow/graph.cpp maxflow/maxflow.cpp)
target_link_libraries(test_graphcut ASTex)

add_executable(test_graphcut_2 test_graphcut_2.cpp graphcut.h graphcut.cpp maxflow/block.h maxflow/graph.h maxflow/graph.cpp maxflow/maxflow.cpp)
target_link_libraries(test_graphcut_2 ASTex)
