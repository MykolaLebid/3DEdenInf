MAIN_CATALOGUE_NAME="all_results/FS_survivor_4_10_num_16052018"    	

SIMULATION_MODE='1' # 1 - for ABC inference; 2 - distance comparison; 3 - statistics 
NUM_OF_POINTS_IN_SEARCH_SPACE="1000"
NUM_OF_REPETITIONS="1"
LEFT_SCALE="0.9"
RIGHT_SCALE="0.9"
MODE_OF_DIST="2" # 1 - iterate all parameters; 2 - iterate the first element (second and third like in 1)  
STATISTICS="averaging"

###compile main sim file
cd sim_prog
 ./compiling.bat
cd ..



# file system creation
  # compiling and run
  case $SIMULATION_MODE in
  '1')
  	cd file_system_creation
  	cd dev	 
   		./compiling.bat   
   		./run.bat $MAIN_CATALOGUE_NAME $SIMULATION_MODE $NUM_OF_POINTS_IN_SEARCH_SPACE $NUM_OF_REPETITIONS $LEFT_SCALE $RIGHT_SCALE $MODE_OF_DIST
  	cd ..
  	cd ..	
  ;;
 '2')
  	cd file_system_creation
  	cd dev	 
  		./compiling.bat   
   		./run.bat $MAIN_CATALOGUE_NAME $SIMULATION_MODE $NUM_OF_POINTS_IN_SEARCH_SPACE $NUM_OF_REPETITIONS $LEFT_SCALE $RIGHT_SCALE $MODE_OF_DIST
  	cd ..
  	cd ..	
  ;;
#  '3')
#  	cd file_system_creation
#  	cd dev	 
#		./compiling.bat   
#  		./run.bat $MAIN_CATALOGUE_NAME $SIMULATION_MODE $STATISTICS
#  	cd ..
#  	cd ..	
#  ;;
  esac




  case $SIMULATION_MODE in
  '1')
# Move created bash files to correspondent directories 
   	chmod +x file_system_creation/dev/bash_script_probes.bat
   	mv file_system_creation/dev/bash_script_probes.bat sim_prog   	
#	chmod +x file_system_creation/dev/bash_script_ABC_inferences.bat
#        mv file_system_creation/dev/bash_script_ABC_inferences.bat multi
  ;;
#  '2')
#	# Move created bash files to correspondent directories 
#   	chmod +x file_system_creation/dev/bash_script_probes.bat
#   	mv file_system_creation/dev/bash_script_probes.bat sim_prog
#
#   	chmod +x file_system_creation/dev/bash_script_compare_dist.bat
#        mv file_system_creation/dev/bash_script_compare_dist.bat sim_prog
#  ;;
#  '3')
#   	chmod +x file_system_creation/dev/bash_script_statistics.bat        
#        mv file_system_creation/dev/bash_script_statistics.bat sim_prog
#  ;;
  esac

	


# Run bash scripts
  case $SIMULATION_MODE in
  '1')
  	cd sim_prog  
   		./bash_script_probes.bat
  	cd ..
  ;;
#  '2')
#  	cd sim_prog  
#  		./bash_script_probes.bat
#  	cd ..
#  ;;
#  '3')
#  ;;
  esac


#  case $SIMULATION_MODE in
#  '1')
#  cd multi 
#    ./bash_script_ABC_inferences.bat 
#  cd ..
#  ;;
#  '2')
#  cd sim_prog
#    ./bash_script_compare_dist.bat
#  cd ..	
#  ;;
#  '3')
#  cd sim_prog
#     ./bash_script_statistics.bat
#  cd ..		
#  ;;
#  esac  

# Delete bash scripts
  case $SIMULATION_MODE in
  '1')
#  cd multi 
#    rm bash_script_ABC_inferences.bat
#  cd ..
  rm sim_prog/bash_script_probes.bat 	
 ;;
#  '2')
#  cd sim_prog
#    rm bash_script_compare_dist.bat 
#  cd ..	
#  rm sim_prog/bash_script_probes.bat 	
#  ;;
#  '3')
#  rm sim_prog/bash_script_statistics.bat
#  ;;
#  esac  


mail -s "Jobs_are_finished" mykola.lebid@bsse.ethz.ch