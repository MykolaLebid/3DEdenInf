cmake_minimum_required(VERSION 3.10) # It will be ok more then 3.0 

# We need C++ comiler and boost
# C we be used as well
project(FileSystemCreation
	VERSION 1.0
	LANGUAGES CXX
)

variable_watch(ALL)
#include_directories(${CMAKE_CURRENT_LIST_DIR}/../../config/source)#TODO think how to change this
#target_include_directories(configTarget PRIVATE ../../config/source)
add_executable(file_system_creation 
	main.cpp
	file_sys_for_parallel_computations.cpp
)

target_include_directories(file_system_creation PRIVATE ../../config/source)
target_link_libraries(file_system_creation configTarget )

#message("INCLUDE_DIRECTORIES:")
#foreach(loopVar IN LISTS ${INCLUDE_DIRECTORIES})
#	message("Dir = ${loopVar}")
#endforeach()	

#list(LENGTH INCLUDE_DIRECTORIES len)
#message("length = ${len}")

#message("top CMAKE_SOURCE_DIR		=${CMAKE_SOURCE_DIR}")
#message("top CMAKE_BINARY_DIR		=${CMAKE_BINARY_DIR}")
#message("top CMAKE_CURRENT_SOURCE_DIR	=${CMAKE_CURRENT_SOURCE_DIR}")
#message("top CMAKE_CURRENT_BINARY_DIR	=${CMAKE_CURRENT_BINARY_DIR}")

#add_subdirectory()
#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.60.00 COMPONENTS filesystem system) 
#
if(Boost_FOUND)
    message("boost is intsalled")
    include_directories(${Boost_INCLUDE_DIRS}) 
#    add_executable(progname file1.cxx file2.cxx) 
    target_link_libraries(file_system_creation  ${Boost_LIBRARIES})
endif()


