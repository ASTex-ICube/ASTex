cmake_minimum_required(VERSION 3.0)
 
project(ASTexTests LANGUAGES CXX)

find_package(ASTex)

add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )


set(exec_list test_exr_io test_super_pixel test_saliency nlutz_tests)

foreach(exec ${exec_list})
	add_executable(${exec} ${exec}.cpp)
        target_link_libraries(${exec} ASTex)
endforeach(exec)

# nice vs solution
set_target_properties(${exec_list} PROPERTIES FOLDER tests)
