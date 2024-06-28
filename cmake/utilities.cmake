
include(CMakeParseArguments)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)

##############################################################################################
#                                 asatex_create_package macro                                 #
# This macro is a helper to create package configuration and version files. These files are  #
# needed when using the find_package command.                                                #
# This macro generate 2 versions of each file : one for the build tree and another for the   #
# install tree.                                                                              #
# Build tree:                                                                                #
# 1.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Targets.cmake             #
# 2.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Config.cmake              #
# 3.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>ConfigVersion.cmake       #
#                                                                                            #
# Install tree:                                                                              #
# 1.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Targets.cmake           #
# 2.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Config.cmake            #
# 3.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>ConfigVersion.cmake     #
#                                                                                            #
# Usage example : find_package(astex);                                                       #
# By convention they have to define the following two variables:                             #
# cmake/<cmake-project-name>_LIBRARIES                                                       #
# cmake/<cmake-project-name>_INCLUDE_DIRS                                                    #
##############################################################################################

cmake_minimum_required(VERSION 3.5)

include(CMakeParseArguments)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)


macro(astex_create_package package_root_dir include_dirs_build_tree include_dirs_install_tree itk_direct)

set(UPPER_NAME "")
string(TOUPPER ${PROJECT_NAME} UPPER_NAME)

set(${UPPER_NAME}_INCLUDE_DIRS ${include_dirs_build_tree})

set(ITK_DIR ${itk_direct})

export(TARGETS ${PROJECT_NAME} FILE "${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake")

configure_package_config_file(
	"${package_root_dir}/${PROJECT_NAME}Config.cmake.in"
	"${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}Config.cmake"
	PATH_VARS ${UPPER_NAME}_INCLUDE_DIRS ITK_DIR
	INSTALL_DESTINATION "${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}"
)

write_basic_package_version_file(
	"${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${ML_VERSION_MAJOR}.${ML_VERSION_MINOR}.${ML_VERSION_PATCH}
	COMPATIBILITY ExactVersion
)


set(${UPPER_NAME}_INCLUDE_DIRS  ${include_dirs_install_tree})

install(TARGETS ${PROJECT_NAME}
	EXPORT ${PROJECT_NAME}Targets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(EXPORT ${PROJECT_NAME}Targets DESTINATION "lib/cmake/${PROJECT_NAME}")

write_basic_package_version_file(
	"${CMAKE_BINARY_DIR}/share/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${CGOGN_VERSION_MAJOR}.${CGOGN_VERSION_MINOR}.${CGOGN_VERSION_PATCH}
	COMPATIBILITY ExactVersion
)

set(CURRENT_LIBRARY "${PROJECT_NAME}")
configure_package_config_file(
	"${package_root_dir}/${PROJECT_NAME}Config.cmake.in"
	"${CMAKE_BINARY_DIR}/share/cmake/${PROJECT_NAME}/${PROJECT_NAME}InstallConfig.cmake"
	PATH_VARS ${UPPER_NAME}_INCLUDE_DIRS ITK_DIR
	INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
)

install(FILES "${CMAKE_BINARY_DIR}/share/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake" DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_BINARY_DIR}/share/cmake/${PROJECT_NAME}/${PROJECT_NAME}InstallConfig.cmake" DESTINATION "lib/cmake/${PROJECT_NAME}" RENAME "${PROJECT_NAME}Config.cmake")

endmacro()

##############################################################################################

