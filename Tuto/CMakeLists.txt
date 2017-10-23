cmake_minimum_required(VERSION 3.0)
 
project(ASTexTutos LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


set( exec_list  tuto_pixel_types tuto_gray tuto_rgb tuto_mask 
 tuto_traverse_iterators tuto_traverse_for_all 
 tuto_color_filters
 tuto_filter1 tuto_filter2 tuto_filter3 tuto_filter4 tuto_filter5 tuto_filter6 tuto_filter7 
 tuto_png_indexed )

foreach(exec ${exec_list})
	add_executable(${exec} ${exec}.cpp)
	target_link_libraries(${exec} ASTex)
endforeach(exec)

# nice vs solution
set_target_properties(${exec_list} PROPERTIES FOLDER tests)

# nice vs solution
set_target_properties(tuto_pixel_types tuto_rgb tuto_gray tuto_filter1 tuto_filter2 tuto_filter3 tuto_filter4 tuto_filter5 tuto_filter6 tuto_filter7 tuto_png_indexed tuto_traverse_iterators tuto_color_filters tuto_mask tuto_traverse_for_all PROPERTIES FOLDER tutorials)
