# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(ccvt_test_Charline)


# Include this package's headers first
#include_directories(BEFORE ./ ./include)
set(CMAKE_CXX_STANDARD 17)




## Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core)

#find_package(ASTex REQUIRED)

#add_subdirectory(glfw)
find_package(OpenGL)


set( SUITESPARSE_LIBRARIES
        -lm -lamd -lcamd -lcolamd -lccolamd -lcholmod -lspqr -ltbb -lmetis -lblas -llapack )

set(SRCS main_ccvt_cgal.cpp
        matrix/sparse_array.cpp
        matrix/sparse_matrix.cpp
        imgui/imgui_impl_opengl3.cpp
        imgui/imgui_impl_glfw.cpp
        imgui/imgui.cpp
        imgui/imgui_demo.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_widgets.cpp
        ccvt_assign.cpp
        ccvt_energy.cpp
        ccvt_histogram.cpp
        ccvt_init.cpp
        ccvt_io.cpp
        ccvt_optimizer.cpp
        ccvt_sites.cpp
        ccvt_application.cpp)

add_executable(ccvt_test_main ${SRCS} )

target_sources(ccvt_test_main PRIVATE main_ccvt_cgal.cpp)

target_link_libraries(ccvt_test_main PRIVATE
        CGAL::CGAL
        ${SUITESPARSE_LIBRARIES}
        glfw
        OpenGL::GL
        OpenGL::GLU
        GLEW)
#target_link_libraries(ccvt_test_main PRIVATE ASTex)

# Link with SuiteSparse
#target_link_libraries( ccvt_test_main ${SUITESPARSE_LIBRARIES} )


cgal_add_compilation_test(ccvt_test_main)

# nice vs solution
set_target_properties(ccvt_test_main PROPERTIES FOLDER Algo)
