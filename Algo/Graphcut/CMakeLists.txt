cmake_minimum_required(VERSION 3.5)
 
project(Graphcut)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


#find_package(ASTex)  ## done outside

add_executable(test_graphcut "")
target_sources(test_graphcut PRIVATE test_graphcut.cpp graphcut.h graphcut.cpp maxflow/block.h maxflow/graph.h maxflow/graph.cpp maxflow/maxflow.cpp)
target_link_libraries(test_graphcut PRIVATE ASTex)

add_executable(test_graphcut_2 "")
target_sources(test_graphcut_2 PRIVATE test_graphcut_2.cpp graphcut.h graphcut.cpp maxflow/block.h maxflow/graph.h maxflow/graph.cpp maxflow/maxflow.cpp)
target_link_libraries(test_graphcut_2 PRIVATE ASTex)

# nice vs solution
set_target_properties(test_graphcut test_graphcut_2 PROPERTIES FOLDER Algo)
