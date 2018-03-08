#include <iostream>
#include "file_sys_for_parallel_computations.h"

int main(int argc, char *argv[])
{
    std::string init_conf_file_name =  "init_conf_file.dat";
    std::string step_lengths_file_name =  "step_lengths.dat";
    std::string step_numbers_file_name =  "step_numbers.dat";
    std::string orders_file_name =  "orders.dat";
    std::string sim_directory_name = argv[1];
	  //creates directories and setting files according to 4 main files (begin)
		StatisticsFS statistics_probe_fs;
		statistics_probe_fs.get(init_conf_file_name,
														step_lengths_file_name,
														step_numbers_file_name,
														orders_file_name);

	  //(end)
		switch (atoi(argv[2])){
			case 1:{
				statistics_probe_fs.createFS(sim_directory_name);
				std::string bash_file_name = "bash_script_probes.bat";
				statistics_probe_fs.createBashEtalonSims(bash_file_name,
																						 sim_directory_name,
																						 init_conf_file_name);
				std::string bash_file_name_ABC = "bash_script_ABC_inferences.bat";
				statistics_probe_fs.createBashForABCInference(bash_file_name_ABC,
																										  sim_directory_name,
																											init_conf_file_name);
			}; break;
			case 2:{
				statistics_probe_fs.createFS(sim_directory_name);
				std::string bash_file_name = "bash_script_probes.bat";
				statistics_probe_fs.createBashEtalonSims(bash_file_name,
																						 sim_directory_name,
																						 init_conf_file_name);
				std::string bash_file_name_compare_dist = "bash_script_compare_dist.bat";
				statistics_probe_fs.createBashCompareDist(bash_file_name_compare_dist,
																									sim_directory_name,
																									init_conf_file_name,
																									atoi(argv[3]),
																									atoi(argv[4]),
																									atof(argv[5]),
																									atof(argv[6]),
																									atoi(argv[7]));
			}; break;
			case 3:{
				std::string mode = argv[3];
				std::string bash_file_name_statistics = "bash_script_statistics.bat";
				statistics_probe_fs.initFSforStatistics(sim_directory_name);
				statistics_probe_fs.createBashStatistics(bash_file_name_statistics,
																								 sim_directory_name,
																								 init_conf_file_name,
																								 mode);
			}; break;
			default:{
				std::cout<< "Template: <catalog directory (will be in upper directory)> <mode: 1 or 2> (for mode 2 <number of points in search space(2*num + 1) > <number of repetitions> <mode: 0 or 1 (0 - 1 iterate one parameter; 1 - iterate all parameters)>)"<< std::endl
				         << "There are 2 possibilities: 1) create file infrastructure for ABC inference" << std::endl
								 << "Example: FS_test 1" <<  std::endl
								 << "2) create file infrastructure for distance comparisons" <<std::endl
								 << "Example: FS_test 2 10 10 0" << std::endl
								 << "3) work with a statistical data" << std::endl
								 << "Example: FS_test 3 averaging" << std::endl;
			}; break;
	  };
    return 0;
}
