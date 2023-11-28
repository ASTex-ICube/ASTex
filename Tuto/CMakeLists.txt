cmake_minimum_required(VERSION 3.0)
 
project(ASTexTutos LANGUAGES CXX)

find_package(ASTex)

set( exec_list 
 tuto_pixel_types tuto_gray tuto_rgb tuto_mask 
 tuto_traverse_iterators tuto_traverse_for_all 
 tuto_color_filters
 #tuto_filter1 tuto_filter2 tuto_filter3 tuto_filter4 tuto_filter5 tuto_filter6 tuto_filter7 
 tuto_png_indexed
 tuto_store
 tuto_mp
 tuto_tnb)

foreach(exec ${exec_list})
	add_executable(${exec} "")
	target_sources(${exec} PRIVATE ${exec}.cpp)
	target_link_libraries(${exec} PRIVATE ASTex)
endforeach(exec)

# nice vs solution
set_target_properties(${exec_list} PROPERTIES FOLDER tutorials)
