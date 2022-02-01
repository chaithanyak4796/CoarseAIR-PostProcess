cmake_minimum_required (VERSION 2.8.11)
project (Statistics)

# Do you want to compile with OpenMP support? (YES/NO) : OpenMP jumbles the log file.
set (Include_OpenMP NO)

set(EXECUTABLE_OUTPUT_PATH ../)
set(CMAKE_SOURCE_DIR ./src/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES

###### Headers #########
    ${CMAKE_SOURCE_DIR}Global.h
    ${CMAKE_SOURCE_DIR}Logger.h
    ${CMAKE_SOURCE_DIR}Main.h
    ${CMAKE_SOURCE_DIR}Input_Class.h
    ${CMAKE_SOURCE_DIR}Statistics.h 
    ${CMAKE_SOURCE_DIR}Reactions.h

###### Functions #######
    ${CMAKE_SOURCE_DIR}Main.cpp
    ${CMAKE_SOURCE_DIR}Find_length.cpp
   
###### Classes #########
    ${CMAKE_SOURCE_DIR}Input_Class.cpp
    ${CMAKE_SOURCE_DIR}Statistics.cpp
    ${CMAKE_SOURCE_DIR}Reactions.cpp
)


add_executable(Statistics ${SOURCES})

if(Include_OpenMP)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	target_link_libraries(CompileTrajs PUBLIC OpenMP::OpenMP_CXX)
endif()	
endif()
