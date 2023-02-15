cmake_minimum_required(VERSION 3.0)

project(micropatterns)

find_package(ASTex REQUIRED)

set (MICRO_PATTERNS_SRC
vec_ops.h
mp_inputs.h
mipmaped_noise.h
blurred_colorMap.h
tiling.h
tiling.cpp
micropatterns.h
)


add_library(${PROJECT_NAME} "")
target_sources(${PROJECT_NAME} PRIVATE ${MICRO_PATTERNS_SRC})
target_link_libraries(${PROJECT_NAME} PUBLIC ASTex)


# nice vs solution
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER Algo)


add_executable(test_mp "")
target_sources(test_mp PRIVATE test_mp.cpp)
target_link_libraries(test_mp PRIVATE ASTex ${PROJECT_NAME})