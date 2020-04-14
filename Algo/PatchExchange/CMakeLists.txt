cmake_minimum_required(VERSION 3.0)
 
project(PatchExchange LANGUAGES CXX)

find_package(ASTex)

add_executable(test_patch_exchange "")

target_sources(test_patch_exchange PRIVATE test_patch_exchange.cpp
	FragmentProcessor.cpp PatchProcessor.cpp PoissonDiskSampler.cpp
	FragmentProcessor.h PatchProcessor.h PoissonDiskSampler.h)

target_link_libraries(test_patch_exchange PRIVATE ASTex)

# nice vs solution
set_target_properties(test_patch_exchange PROPERTIES FOLDER tests)
