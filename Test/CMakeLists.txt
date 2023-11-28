cmake_minimum_required(VERSION 3.0)
 
project(ASTexTests LANGUAGES CXX)

find_package(ASTex)
#find_package(Qt5 REQUIRED COMPONENTS Core Widgets PrintSupport)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )

set(exec_list test_exr_io test_super_pixel test_saliency test_contentExchange test_sparse
    test_matchHistogram test_RPN test_terrainAmplification test_gradient test_multiTextureContentExchange
    test_heightToTexture test_PCTSResynthesis test_texton test_contentExchangeAmplification test_csn 
	test_spectrum test_cyclostationaryPSD
	test_simple_tiling_n_blending
	test_mixmax)


#find_path(HeightFieldLib_path REQUIRED)
#find_library(HeightFieldLib REQUIRED)
#find_path(CoreLib_path REQUIRED)
#find_library(CoreLib REQUIRED)
#find_path(MayaLib_path REQUIRED)
#find_library(MayaLib REQUIRED)



#include_directories(${HeightFieldLib_path} ${CoreLib_path} ${MayaLib_path})
#set(LIBS ${LIBS} ${HeightFieldLib} ${CoreLib} ${MayaLib} Qt5::Core Qt5::Widgets Qt5::PrintSupport)

foreach(exec ${exec_list})
	add_executable(${exec} "")
	target_sources(${exec} PRIVATE ${exec}.cpp)
	target_link_libraries(${exec} PRIVATE ASTex ${LIBS})
endforeach(exec)


# nice vs solution
set_target_properties(${exec_list} PROPERTIES FOLDER tests)
