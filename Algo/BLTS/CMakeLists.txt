cmake_minimum_required(VERSION 3.0)
 
project(LRPN LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


file(GLOB HV_HEADERS hv/*.h)

add_executable(blt_synthesis "")
target_sources(blt_synthesis PRIVATE
	blt_synthesis.cpp
    noise_filter.cpp
    noise_kmean.cpp
    noise_synthesis_correlated_color.cpp
    almost_wang_tiles.cpp biscalenoisepatchexg.cpp
    hv_astex_io.h ${HV_HEADERS})

target_link_libraries(blt_synthesis PRIVATE ASTex)

add_executable(test_noise_filter "")
target_sources(test_noise_filter PRIVATE test_noise_filter.cpp noise_filter.cpp)
target_link_libraries(test_noise_filter ASTex)

add_executable(test_noise_kmean "")
target_sources(test_noise_kmean PRIVATE test_noise_kmean.cpp noise_kmean.cpp)
target_link_libraries(test_noise_kmean PRIVATE ASTex)

add_executable(test_noise_synthesis_correlated_color "")
target_sources(test_noise_synthesis_correlated_color PRIVATE
	test_noise_synthesis_correlated_color.cpp noise_synthesis_correlated_color.cpp)
target_link_libraries(test_noise_synthesis_correlated_color PRIVATE ASTex)


# nice vs solution
set_target_properties(
test_noise_kmean
test_noise_synthesis_correlated_color
test_noise_filter
blt_synthesis 
PROPERTIES FOLDER LRPN)
