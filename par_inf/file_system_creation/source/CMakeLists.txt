cmake_minimum_required(VERSION 3.3) # It will be ok more then 3.0 

# We need C++ comiler and boost
# C we be used as well
project(FileSystemCreation
	VERSION 1.0
	LANGUAGES CXX
)

add_executable(file_system_creation 
	main.cpp
	file_sys_for_parallel_computations.cpp
)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.60.00 COMPONENTS *boost libraries here*) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    add_executable(progname file1.cxx file2.cxx) 
    target_link_libraries(progname ${Boost_LIBRARIES})
endif()


