cmake_minimum_required(VERSION 3.5)

project(ASTex CXX C)

#if (MSVC)
#find_package(GTest REQUIRED)
#endif()

set(ASTEX_VERSION_MAJOR 0)
set(ASTEX_VERSION_MINOR 9)
set(ASTEX_VERSION_PATCH 0)

#cmake_policy(SET CMP0167 NEW)


set(ENV{PATH}=${CMAKE_PREFIX_PATH};$ENV{PATH})

#set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH  ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
include ( ${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

IF(WIN32)
STRING(REGEX REPLACE "\\\\" "/" TEMPO_P $ENV{userprofile}/Desktop/ASTex_data)
	set(ASTEX_TEMPO_PATH ${TEMPO_P} CACHE STRING "tempo directory (tutos img) storage" )
ELSE()
	set(ASTEX_TEMPO_PATH /tmp/ASTex_data CACHE STRING "tempo directory (tutos img) storage" )
ENDIF()
file (GLOB DATA_CONTENTS "${CMAKE_SOURCE_DIR}/Data/*")
MESSAGE(WARNING "ASTEX TEMP PATH = " ${ASTEX_TEMPO_PATH})
IF(NOT EXISTS ${ASTEX_TEMPO_PATH})
    MESSAGE(WARNING "NO ASTEX TEMP PATH = ")
    MESSAGE(WARNING ${DATA_CONTENTS})
    file(COPY ${DATA_CONTENTS} DESTINATION "${ASTEX_TEMPO_PATH}")
ENDIF()


option(ASTEX_BUILD_Bench "compile benchmark" ON )
option(ASTEX_BUILD_Test "compile tests" ON )
option(ASTEX_BUILD_Tuto "compile tutorials" ON )
option(BUILD_SHARED_LIBS "build dynamic lib" OFF )
#option(ASTEX_USE_CPP14 "set this to on if VXL say that you are using a C++ standard version older than ..." OFF )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(ASTEX_SOURCE_DIR ${CMAKE_SOURCE_DIR})


IF(NOT WIN32)
# deduce Release/Debug from build directory name (except in Windows)
  IF (NOT DEFINED CMAKE_BUILD_TYPE)
    IF (${CMAKE_CURRENT_BINARY_DIR} MATCHES "(.*)Debug|(.*)debug|(.*)DEBUG")
      SET(CMAKE_BUILD_TYPE Debug)
    ELSE()
      SET(CMAKE_BUILD_TYPE Release)
    ENDIF()
  ENDIF()
ENDIF()


add_definitions(-DTEMPO_PATH=std::string\(\"${ASTEX_TEMPO_PATH}/\"\) )

add_subdirectory( ${ASTEX_SOURCE_DIR}/ASTex )

if (${ASTEX_BUILD_Bench})
	add_subdirectory( Bench )
endif()

if (${ASTEX_BUILD_Test})
        add_subdirectory( Test )
endif()

if (${ASTEX_BUILD_Tuto})
        add_subdirectory( Tuto )
endif()

if (MSVC)
	add_subdirectory( deploy_dll )
endif()


set ( ALGO_DIR ${ASTEX_SOURCE_DIR}/Algo)
file(GLOB algos RELATIVE ${ALGO_DIR} ${ALGO_DIR}/*)

#find_package(ASTex)  ## done outside
foreach(child ${algos})
  if(IS_DIRECTORY  ${ALGO_DIR}/${child})
    if(EXISTS "${ALGO_DIR}/${child}/CMakeLists.txt")
      option(ASTEX_ALGO_${child} "compile $(child}" ON )
    endif()
  endif()
endforeach()

foreach(child ${algos})
  if(IS_DIRECTORY  ${ALGO_DIR}/${child})
    if(EXISTS "${ALGO_DIR}/${child}/CMakeLists.txt")
      if (${ASTEX_ALGO_${child}})
        add_subdirectory( Algo/${child} )
      endif()
    endif()
  endif()
endforeach()

set(EXCLUDED_DIRS "ASTex|Tuto|Test|Bench|Algo|cmake|Data|Eigen|scrip_windows|deploy_dll")


file(GLOB children RELATIVE  ${ASTEX_SOURCE_DIR} ${ASTEX_SOURCE_DIR}/*)
foreach(child ${children})
  if(IS_DIRECTORY  ${ASTEX_SOURCE_DIR}/${child})
    if(EXISTS "${ASTEX_SOURCE_DIR}/${child}/CMakeLists.txt")
      if (NOT("${child}"  MATCHES "^(${EXCLUDED_DIRS})$"))
        option(ASTEX_PERSO_${child} "compile $(child}" OFF )
      endif()
    endif()
  endif()
endforeach()

foreach(child ${children})
  if(IS_DIRECTORY  ${ASTEX_SOURCE_DIR}/${child})
    if(EXISTS "${ASTEX_SOURCE_DIR}/${child}/CMakeLists.txt")
      if (NOT("${child}"  MATCHES "^(${EXCLUDED_DIRS})$"))
        if (${ASTEX_PERSO_${child}})
          add_subdirectory( ${child} )
        endif()
      endif()
    endif()
  endif()
endforeach()


