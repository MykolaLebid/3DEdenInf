###compile main sim file

cd sim_prog
# ./compiling.bat
cd ..

cd multi
# ./compiling.bat
cd ..

# file system creation
  # compiling and run
  cd file_system_creation
  cd dev	 
   ./compiling.bat   
   ./run.bat    	
  cd ..
  cd ..	

# Move created bash files to correspondent directories 
  chmod +x file_system_creation/dev/bash_script_probes_creation.bat
  chmod +x file_system_creation/dev/bash_script_ABC_inferences.bat

  mv file_system_creation/dev/bash_script_probes_creation.bat sim_prog
  mv file_system_creation/dev/bash_script_ABC_inferences.bat multi

# Run bush scripts
  cd sim_prog
    ./bash_script_probes_creation.bat 
  cd ..
  
  cd multi 
    ./bash_script_ABC_inferences.bat 
  cd ..

# correspondent directories
#  rm sim_prog/bash_script_probes_creation.bat
#  rm multi/bash_script_ABC_inferences.bat
