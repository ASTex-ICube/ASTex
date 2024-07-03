cmake_minimum_required(VERSION 3.1...3.23)
project(ccvt_cgal_chgrenier)

## Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core)

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)

if (WIN32)
	find_package(GLEW CONFIG REQUIRED)
	find_package(metis CONFIG REQUIRED)
	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)
	find_package(TBB CONFIG REQUIRED)
	find_package(suitesparse CONFIG REQUIRED)

	set( SUITESPARSE_LIBRARIES
        	${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
		SuiteSparse::amd
		SuiteSparse::camd
		SuiteSparse::colamd
		SuiteSparse::ccolamd
	    SuiteSparse::cholmod
		SuiteSparse::spqr 
		TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy )
	set (GLEW_LIBRARIES GLEW::GLEW)
else ()

	set( SUITESPARSE_LIBRARIES
        	${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
        	amd camd colamd ccolamd cholmod spqr 
		tbb )
	set (GLEW_LIBRARIES GLEW)
endif()
	
set(SRCS main_ccvt_cgal.cpp
        # matrix
        matrix/sparse_array.cpp
        matrix/sparse_matrix.cpp
        #imgui
        imgui/imgui_impl_opengl3.cpp
        imgui/imgui_impl_glfw.cpp
        imgui/imgui.cpp
        imgui/imgui_demo.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_widgets.cpp
        imgui/imgui_tables.cpp
        #implot
        imgui/implot.cpp
        imgui/implot_items.cpp
        imgui/implot_demo.cpp
        # ccvt
        ccvt_assign.cpp
        ccvt_energy.cpp
        ccvt_histogram.cpp
        ccvt_init.cpp
        ccvt_io.cpp
        ccvt_optimizer.cpp
        ccvt_sites.cpp
        ccvt_application.cpp)

add_executable(ccvt_cgal_chgrenier )

target_sources(ccvt_cgal_chgrenier PRIVATE ${SRCS} main_ccvt_cgal.cpp)

target_include_directories(ccvt_cgal_chgrenier PRIVATE 
				/usr/include/suitesparse
				/usr/local/include/suitesparse)

target_link_directories(ccvt_cgal_chgrenier PRIVATE /usr/local/lib/)

target_link_libraries(ccvt_cgal_chgrenier PRIVATE
        CGAL::CGAL
        ${SUITESPARSE_LIBRARIES}
        glfw
        OpenGL::GL
        OpenGL::GLU
	${GLEW_LIBRARIES} )



target_compile_features(ccvt_cgal_chgrenier PUBLIC cxx_std_17)

cgal_add_compilation_test(ccvt_cgal_chgrenier)

# nice vs solution
set_target_properties(ccvt_cgal_chgrenier PROPERTIES FOLDER Algo)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION "${ASTEX_TEMPO_PATH}")


#install(TARGETS ccvt_cgal_chgrenier
#        COMPONENT ink_edit_comp
#        DESTINATION ink_edit
#)
#install(IMPORTED_RUNTIME_ARTIFACTS ccvt_cgal_chgrenier
#        RUNTIME_DEPENDENCY_SET ink_edit_dep
#        RUNTIME
#        COMPONENT ink_edit_comp
#        DESTINATION ink_edit/lib
#)
#install(RUNTIME_DEPENDENCY_SET ink_edit_dep
#        COMPONENT ink_edit_comp
#        DESTINATION ink_edit/lib
#)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shaders
#        COMPONENT ink_edit_comp
#        DESTINATION ink_edit
#)
#install(CODE [[
#        file(GET_RUNTIME_DEPENDENCIES
#            RESOLVED_DEPENDENCIES_VAR ink_edit_dep_var
#            EXECUTABLES $<TARGET_FILE:ccvt_cgal_chgrenier>
#        )
#        MESSAGE(STATUS "blabla ${ink_edit_dep_var}")
#        file(INSTALL ${ink_edit_dep_var}
#            DESTINATION ink_edit/lib
#        )
#        ]]
#        COMPONENT ink_edit_comp
#)
#install(FILES ${ink_edit_dep_var}
#        COMPONENT ink_edit_comp
#        DESTINATION ink_edit/lib
#)

set_target_properties(ccvt_cgal_chgrenier PROPERTIES INSTALL_RPATH $ORIGIN/lib)
