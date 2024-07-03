cmake_minimum_required(VERSION 3.5)

project(ProceduralNoiseFiltering)

find_package(ASTex REQUIRED)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )

#add_executable(pnf_texture "")
#target_sources(pnf_texture PRIVATE pnf_texture.cpp
#    pnf.h
#    texture_noise.h
#    color_map.h
#    gaussian_transfer.h
#    histogram.h)
#target_link_libraries(pnf_texture ASTex)

#add_executable(pnf "")
#target_sources(pnf PRIVATE pnf.cpp
#    pnf.h
#    noise.h
#    color_map.h
#    gaussian_transfer.h
#    gaussian_transfer.cpp)
#target_link_libraries(pnf ASTex)


add_executable(generator "")
target_sources(generator PRIVATE generator.cpp
	histogram.h
	gaussian_transfer.h
	color_map.h
	pnf.h
	texture_noise.h
	noise.h
	utils.h)
target_link_libraries(generator ASTex)
set_target_properties(generator PROPERTIES FOLDER Algo/ProceduralNoiseFiltering)
