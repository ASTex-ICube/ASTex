cmake_minimum_required(VERSION 3.0)
 
project(PatchExchange LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


set (HEADERS FragmentProcessor.h PatchProcessor.h PoissonDiskSampler.h)
set (SOURCES FragmentProcessor.cpp PatchProcessor.cpp PoissonDiskSampler.cpp)


add_executable(test_patch_exchange test_patch_exchange.cpp  ${SOURCES} ${HEADERS})
target_link_libraries(test_patch_exchange ASTex)

