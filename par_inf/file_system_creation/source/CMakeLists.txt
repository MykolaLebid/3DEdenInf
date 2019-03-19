cmake_minimum_required(VERSION 3.10) # It will be ok more then 3.0 

# We need C++ comiler and boost
# C we be used as well
project(FileSystemCreation
	VERSION 1.0
	LANGUAGES CXX
)

#variable_watch(ALL)
#include_directories(${CMAKE_CURRENT_LIST_DIR}/../../config/source)#TODO think how to change this
#target_include_directories(configTarget PRIVATE ../../config/source)
add_executable(file_system_creation 
	main.cpp
	file_sys_for_parallel_computations.cpp
)

target_include_directories(file_system_creation PRIVATE ../../config/source)
target_link_libraries(file_system_creation configTarget)

#message("TARGET_:")
#foreach(loopVar IN LISTS ${TARGET_LIBRARIES})
#	message("Dir = ${loopVar}")
#endforeach()	

#list(LENGTH TARGET_OBJECTS len)
#message("length = ${len}")

#get_property(properyValue TARGET PROPERTY configTarget BRIEF_DOCS)
#message("${properyValue}")

#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS configTarget PROPERTIES TYPE INCLUDE_DIRECTORIES LINK_LIBRARIES  RESOURCE PUBLIC_HEADER PRIVATE_HEADER
#	SOURCE_DIR VERSION NAME)
#cmake_print_properties(TARGETS file_system_creation PROPERTIES TYPE COMPILE_OPTIONS HEADER_FILE_ONLY INCLUDE_DIRECTORIES COMPILE_FLAGS COMPILE_OPTIONS LINK_DEPENDS LINKER_LANGUAGE LINK_FLAGS LINK_DIRECTORIES LINK_OPTIONS LINK_LIBRARIES ABSTRACT COMPILE_DEFINITIONS LANGUAGE)

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
    #include_directories(${Boost_INCLUDE_DIRS}) 
    #message("Boost_INCLUDE_DIRS")
    #foreach(loopVar IN LISTS ${Boost_INCLUDE_DIRS})
    #	message("Dir = ${loopVar}")
    #endforeach()	



    target_include_directories(file_system_creation PUBLIC ${Boost_INCLUDE_DIRS})
    #    add_executable(progname file1.cxx file2.cxx) 
    target_link_libraries(file_system_creation  ${Boost_LIBRARIES})
endif()

#cmake_print_properties(TARGETS file_system_creation PROPERTIES TYPE COMPILE_OPTIONS HEADER_FILE_ONLY INCLUDE_DIRECTORIES COMPILE_FLAGS COMPILE_OPTIONS LINK_DEPENDS LINKER_LANGUAGE LINK_FLAGS LINK_DIRECTORIES LINK_OPTIONS LINK_LIBRARIES ABSTRACT COMPILE_DEFINITIONS LANGUAGE C_STANDARD C_EXTENSIONS C_STANDARD_REQUIRED HAS_CXX OUTPUT_NAME)

#add_custom_command(TARGET file_system_creation POST_BUILD COMMAND  )
