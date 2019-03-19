add_subdirectory(source) 
#target_include_directories(configTarget PRIVATE source)
#include_directories(	${CMAKE_CURRENT_LIST_DIR}/source)
#set(copyDir "${CMAKE_CURRENT_BINARY_DIR}/settings")
#set(copyDir "${CMAKE_CURRENT_BINARY_DIR}/settings")
#add_custom_target(copySettingFiles 
#	          COMMAND ${CMAKE_COMMAND} -E echo "FACK UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#	          COMMAND ${CMAKE_COMMAND} -E copy_directory "${copyDir}")
#execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${copyDir}")	  
#execute_process(
#	#COMMAND ${CMAKE_COMMAND} -E copy_directory "${copyDir}"
#	COMMAND ${CMAKE_COMMAND} -E echo "FACK UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#	)


file(COPY settings/ DESTINATION settings FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE)
#add_subdirectory(settings)

