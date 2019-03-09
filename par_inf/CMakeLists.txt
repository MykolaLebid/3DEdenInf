cmake_minimum_required(VERSION 3.2) # It will be ok more then 3.0 

# We need C++ comiler and boost
# C we be used as well
project(FileSystemCreation
	VERSION 1.0
	LANGUAGES C CXX
)

add_executable(file_system_creation 
	main.cpp
	file_system_for_parallel_computaion.cpp
)




