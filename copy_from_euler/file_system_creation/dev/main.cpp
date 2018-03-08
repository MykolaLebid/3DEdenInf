#include <iostream>
#include "file_sys_for_parallel_computations.h"
//using namespace std;
//#include <ctime>

int main(int argc, char *argv[])
{
		//std::clock_t start; // Del
		//double duration; // Del
		//start = std::clock(); // Del
		//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; // Del
	  //std::cout<<"time for  (1): "<< duration <<'\n'; // Del

//		ParameterStructClass parameters;
    StatisticsFS statistics_probe_fs;

    //std::string etalon_probe_settings = "etalon_probe_settings";

    std::string init_conf_file_name =  "init_conf_file.dat";
    std::string step_lengths_file_name =  "step_lengths.dat";
    std::string step_numbers_file_name =  "step_numbers.dat";
    std::string orders_file_name =  "orders.dat";

    std::string sim_directory_name = argv[1];

//    parameters.GetFromFile(init_conf_file_name);
//    parameters.PrintToConsole();
		statistics_probe_fs.Get(init_conf_file_name,
														step_lengths_file_name,
														step_numbers_file_name,
														orders_file_name);
		std::string bash_file_name = "bash_script_probes_creation.bat";
//
	  statistics_probe_fs.CreateFS(sim_directory_name);
		statistics_probe_fs.CreateBashForIndependenSim(bash_file_name,
																						       sim_directory_name,
																						       init_conf_file_name);

		std::string bash_file_name_ABC = "bash_script_ABC_inferences.bat";
		statistics_probe_fs.CreateBashForABCInference(bash_file_name_ABC,
																									sim_directory_name,
																									init_conf_file_name);


//   statistics_fs.GetBashForDirInAllCreation(bash_file_name_for_add_dirs,
//																				 sim_derectory_name, additional_dirs);
//		std::vector <std::string> path_vec;
//		path_vec = statistics_probe_fs.GetPathVector(sim_directory_name);

//////		for (const auto & str : path_vec) // access by const reference
//////        std::cout << str << ' ';
//////				std::cout << '\n';
////

////
////		bash_file_name = "bash_script_inference.bat";
////
////		for (auto & str : path_vec) {// access by const reference
////			statistics_inference_fs.CreateFS(str);
////			//statistics_inference_fs.CreateBashForIndependenSim(bash_file_name,
////			//																								   str,
////			//																								   exe_file_name);
////		};
//
    return 0;
}
