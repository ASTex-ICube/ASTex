cmake_minimum_required(VERSION 3.5)

project(ASTex LANGUAGES CXX C)

find_package(PNG REQUIRED)
find_package(ITK CONFIG REQUIRED)
include(${ITK_USE_FILE})
find_package(OpenMP)
find_package(OpenMP)

find_package(Eigen3 REQUIRED)

find_package(Imath CONFIG)
if(NOT TARGET Imath::Imath)
    find_package(IlmBase)
endif()
find_package(OpenEXR CONFIG REQUIRED)

if(__APPLE__)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	set (APPLE_LINK "${CORE_FOUNDATION} ${CARBON}")
endif()


add_library( ${PROJECT_NAME} )


target_sources (${PROJECT_NAME} PRIVATE 
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
	fourier.h
	fourier.hpp
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
	store.h
	slic.h
	region_growing/connected_components.h
	region_growing/neighborhood.h
	region_growing/region_growing.h
	saliency/contrast_based_saliency.h
        saliency/permutohedral_lattice.h

        Sparse/dictionary.h
	Sparse/dictionaryProcessor.h

        mipmap.h
        histogram.h
        rpn_utils.h
        texton_io.h
        imageviewer.h

	../Algo/PatchExchange/FragmentProcessor.h
	../Algo/PatchExchange/PatchProcessor.h
	../Algo/PatchExchange/PoissonDiskSampler.h

        ContentExchange/atlas.h
        ContentExchange/content.h
        ContentExchange/patch.h
        ContentExchange/patchProcessor.h
	ContentExchange/multiTextureProcessor.h
        Stamping/sampler.h
        Stamping/stamp.h
        Stamping/stamper.h
        PCTS/pcts.h
	#filters.h
	colorMatching.h
	RegionFinder/exploreTextures.h
	RegionFinder/comparePixels.h
	RegionFinder/regionFinder.h
	texture_pool.h

	CSN/csn_texture.h
	CSN/tpsd.h

	TilingBlending/simple_tnb.h
)

target_sources (${PROJECT_NAME} PRIVATE
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
	thread_pool.cpp

	../Algo/PatchExchange/FragmentProcessor.cpp
	../Algo/PatchExchange/PatchProcessor.cpp
	../Algo/PatchExchange/PoissonDiskSampler.cpp

        texton_io.cpp
        Stamping/sampler.cpp
        rpn_utils.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC ${ITK_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS}  ${FFTW_INCLUDES} ${OPENEXR_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${ASTEX_SOURCE_DIR}>)



set ( ASTEX_LIBRARIES ${ITK_LIBRARIES} ${FFTW_LIBRARIES} ${PNG_LIBRARIES} ${APPLE_LINK})

# if (NOT(MSVC OR APPLE))
#   target_compile_options(${PROJECT_NAME} PUBLIC -pthread)
# endif()

if(OpenMP_CXX_FOUND)
	target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -D_USE_MATH_DEFINES -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -D_SECURE_SCL=0 -D_SECURE_SCL_THROWS=0)
	if(${BUILD_SHARED_LIBS})
		target_compile_definitions(${PROJECT_NAME} PRIVATE "-DASTEX_DLL_EXPORT")
	else()
		target_compile_definitions(${PROJECT_NAME} PUBLIC "-DBUILD_STATIC")
	endif()
else()
	target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wno-unknown-pragmas -Wno-missing-braces -fPIC )
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC "${PNG_DEFINITION}")
  
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

target_link_libraries(${PROJECT_NAME} PUBLIC ${ASTEX_LIBRARIES})

# use itkeigen if exist (itk5.2)
list(GET ITK_INCLUDE_DIRS 0 ITK_PATH)
if(EXISTS ${ITK_PATH}/itkeigen)
	target_include_directories(${PROJECT_NAME} BEFORE PUBLIC ${ITK_PATH}/itkeigen/)
else()
	target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
endif()

# link with good libs function of openEXR version
target_link_libraries (${PROJECT_NAME} PUBLIC 
        # For OpenEXR/Imath 3.x:
          $<$<TARGET_EXISTS:OpenEXR::OpenEXR>:OpenEXR::OpenEXR>
          $<$<TARGET_EXISTS:Imath::Imath>:Imath::Imath>
          $<$<TARGET_EXISTS:Imath::Half>:Imath::Half>
        # For OpenEXR 2.4/2.5:
          $<$<TARGET_EXISTS:OpenEXR::IlmImf>:OpenEXR::IlmImf>
          $<$<TARGET_EXISTS:IlmBase::Imath>:IlmBase::Imath>
          $<$<TARGET_EXISTS:IlmBase::Half>:IlmBase::Half>
          $<$<TARGET_EXISTS:IlmBase::IlmThread>:IlmBase::IlmThread>
          $<$<TARGET_EXISTS:IlmBase::Iex>:IlmBase::Iex>
    )


file(GLOB HEADERS "." "*.h" "*.hpp")
install(FILES ${HEADERS} DESTINATION include/ASTex/)
install(FILES "${ASTEX_SOURCE_DIR}/cmake/FindOpenEXR.cmake" "${ASTEX_SOURCE_DIR}/cmake/FindFFTW.cmake" DESTINATION lib/cmake/ASTex/)

astex_create_package("${CMAKE_CURRENT_SOURCE_DIR}" "${ASTEX_SOURCE_DIR}" "include" ${ITK_DIR})
