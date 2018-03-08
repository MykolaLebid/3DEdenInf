#include <iostream>
#include "file_sys_for_parallel_computations.h"
using namespace std;
//#include <ctime>


int main(int argc, char *argv[])
{
		//std::clock_t start; // Del
		//double duration; // Del
		//start = std::clock(); // Del
		//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; // Del
	  //std::cout<<"time for  (1): "<< duration <<'\n'; // Del

		Parameters parameters;
    StatisticsFS statistics_fs;
    std::string init_conf_file_name =  "init_conf_file.dat";
    std::string step_lengths_file_name = "step_lengths.dat";
    std::string step_numbers_file_name = "step_numbers.dat";
    std::string orders_file_name =  "orders.dat";
    std::string sim_directory_name = argv[1];
    std::string exe_file_name = argv[2];


    //parameters.GetFromFile(init_conf_file_name);
    //parameters.PrintToConsole();

		statistics_fs.Get( init_conf_file_name, step_lengths_file_name,
											 step_numbers_file_name, orders_file_name);
		std::string bash_file_name = "bash_script_probes_creation.bat";
		std::string bash_file_name_for_add_dirs = "bash_script_for_add_dirs.bat";
		//vector<std::string> additional_dirs;
		//additional_dirs.push_back("BirthRateGrids");
		//additional_dirs.push_back("DriverGrids");
		//additional_dirs.push_back("PassengersGrids");
		//additional_dirs.push_back("PopsGrids");

		statistics_fs.CreateFS(sim_directory_name);
		statistics_fs.CreateBashForIndependenSim(bash_file_name,
																						 sim_directory_name,
																						 exe_file_name);
		//statistics_fs.GetBashForDirInAllCreation(bash_file_name_for_add_dirs,
		//																				 sim_derectory_name, additional_dirs);

    return 0;
}
