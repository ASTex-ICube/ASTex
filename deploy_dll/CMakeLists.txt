cmake_minimum_required(VERSION 3.0)
 
project(deploy_astex_dll)

string(REPLACE "/cl.exe" "/" MSVC_PATH ${CMAKE_CXX_COMPILER})
find_program(DUMPBIN dumpbin.exe PATHS ${MSVC_PATH})

set(DUMPBIN ${DUMPBIN} CACHE FILEPATH "dumpbin exec" )
set(TEMP_DIR $ENV{userprofile}/AppData/Local/Temp CACHE PATH "temp directory" )

add_executable(${PROJECT_NAME} deploy_astex_dll.cpp)

target_compile_definitions(${PROJECT_NAME} PRIVATE -DTEMP_DIR=std::string\(\"${TEMPO_PATH}\"\) )
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDUMPBIN=std::string\(\"${DUMPBIN}\"\) )

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)