cmake_minimum_required(VERSION 3.0)

project(ASTex LANGUAGES CXX)

find_package(OpenEXR REQUIRED)
find_package(PNG REQUIRED)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(OpenMP QUIET)


set ( ASTEX_LIBRARIES ${ITK_LIBRARIES} ${FFTW_LIBRARIES} ${PNG_LIBRARIES} ${OPENEXR_LIBRARIES} )
if (NOT(MSVC OR APPLE))
	set ( ASTEX_LIBRARIES ${ASTEX_LIBRARIES} -pthread )
endif()

set (LIB_INCLUDE
	internal.h
	image_common.h
	image_gray.h
	image_rgb.h
	image_rgba.h
	image_indexed.h
	image_spectral.h
	conversion.h
	png_indexed_io.h
	colorspace_filters.h
	mask.h
	easy_io.h
	fourier.h fourier.hpp
	region_traversor.h
	local_spectrum.h
	utils.h
	special_io.h
	exr_io.h
	image_collector.h
	bilateral_filters.h
	distances_maps.h
	pca.h
	dll.h
	thread_pool.h
	slic.h
	region_growing/connected_components.h
	region_growing/neighborhood.h
	region_growing/region_growing.h
	saliency/contrast_based_saliency.h
        saliency/permutohedral_lattice.h

#        Sparse/atom.h
#        Sparse/dictionary.h

        mipmap.h
        histogram.h
        rpn_utils.h
        texton_io.h
        imageviewer.h
        ContentExchange/atlas.h
        ContentExchange/content.h
        ContentExchange/patch.h
        ContentExchange/patchProcessor.h
        ContentExchange/old_PatchProcessor.h
        ContentExchange/old_FragmentProcessor.h
        ContentExchange/old_PoissonDiskSampler.h
        ContentExchange/benchmarker.h
        Stamping/sampler.h
        Stamping/stamp.h
        Stamping/stamper.h
        PCTS/pcts.h
        filters.h
)

set (LIB_SRC
	png_indexed_io.cpp
	mask.cpp
	easy_io.cpp
	fourier.cpp
	local_spectrum.cpp
	file_paths.cpp
	special_io.cpp
	exr_io.cpp
	image_collector.cpp
	bilateral_filters.cpp
	distances_maps.cpp
	pca.cpp
	thread_pool.cpp

#        Sparse/atom.cpp
#        Sparse/dictionary.cpp

        ContentExchange/old_PatchProcessor.cpp
        ContentExchange/old_FragmentProcessor.cpp
        ContentExchange/old_PoissonDiskSampler.cpp
        ContentExchange/benchmarker.cpp

        texton_io.cpp
        Stamping/sampler.cpp
        rpn_utils.cpp
)

add_library( ${PROJECT_NAME} ${LIB_INCLUDE} ${LIB_SRC})

if(MSVC)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -D_USE_MATH_DEFINES -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -D_SECURE_SCL=0 -D_SECURE_SCL_THROWS=0)
  if(${BUILD_SHARED_LIBS})
    target_compile_definitions(${PROJECT_NAME} PRIVATE "-DASTEX_DLL_EXPORT")
  else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC "-DBUILD_STATIC")
  endif()
else()
  if (ITK_VERSION VERSION_LESS 4.10.0)
    target_compile_options(${PROJECT_NAME} PUBLIC -fpermissive )
    message (INFO " Old ITK version, add -fpermissive compil option")
  endif()

  if (${OpenMP_CXX_FLAGS})
      target_compile_options(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
  endif()
#  target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wno-unknown-pragmas -Wno-missing-braces -Wno-deprecated-register -fPIC )
target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wno-unknown-pragmas -Wno-missing-braces -fPIC )
  if (${ASTEX_USE_CPP14})
      target_compile_options(${PROJECT_NAME} PUBLIC -std=c++14 )
  else()
      target_compile_options(${PROJECT_NAME} PUBLIC -std=c++11 )
  endif()
endif()


target_compile_definitions(${PROJECT_NAME} PUBLIC "${PNG_DEFINITION}")


#set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")


target_include_directories( ${PROJECT_NAME}  BEFORE PUBLIC ${ITK_INCLUDE_DIRS} )

target_include_directories(${PROJECT_NAME} PUBLIC
$<INSTALL_INTERFACE:include>
${PNG_INCLUDE_DIRS} ${OPENEXR_INCLUDE_PATHS} ${FFTW_INCLUDES})

# put local ASTex separate and before to ensure use of local before installed
target_include_directories(${PROJECT_NAME} BEFORE PUBLIC
$<BUILD_INTERFACE:${ASTEX_SOURCE_DIR}>)

target_link_libraries(${PROJECT_NAME} PUBLIC ${ASTEX_LIBRARIES})


file(GLOB HEADERS "." "*.h" "*.hpp")
install(FILES ${HEADERS} DESTINATION include/ASTex/)
install(FILES "${ASTEX_SOURCE_DIR}/cmake/FindOpenEXR.cmake" "${ASTEX_SOURCE_DIR}/cmake/FindFFTW.cmake" DESTINATION lib/cmake/ASTex/)

astex_create_package("${CMAKE_CURRENT_SOURCE_DIR}" "${ASTEX_SOURCE_DIR}" "include" ${ITK_DIR})
