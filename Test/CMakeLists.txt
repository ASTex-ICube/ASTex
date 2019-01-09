cmake_minimum_required(VERSION 3.0)
 
project(ASTexTests LANGUAGES CXX)

find_package(ASTex)

set(exec_list
	test_exr_io
	test_super_pixel
	test_saliency
	nlutz_tests
  test_contentExchange 
  test_pcts
	)

foreach(exec ${exec_list})
	add_executable(${exec} "")
	target_sources(${exec} PRIVATE ${exec}.cpp)
	target_link_libraries(${exec} PRIVATE ASTex)
endforeach(exec)

# nice vs solution
set_target_properties(${exec_list} PROPERTIES FOLDER tests)
