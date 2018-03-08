#ifndef FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED
#define FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED

#include <conf_structures.h>
#include <boost/lexical_cast.hpp>

#define BOOST_SYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <string>
#include <iostream>
//void SetParameters(Parameters & parameters, std::string file_name);
//Parameters & GetParameters(std::string file_name);
// PrintToConsole(Parameters & parameters);

typedef std::tuple<std::string /*full path in ptree*/,
									 std::string /*short path*/,
			 	           float /*step length*/,
			 	           int /*step number*/,
				           int /*order in tree hierarchy*/> multi_type;
//TODO
class StatisticsFS{
public:
  Parameters initial_config; // Structure with all initial configurations.
														 //	This structure will be as a template for all
														 // config files.
	std::string initial_config_f_n; // f_n is file_name
  Parameters step_lengths;   // Step length for every numerical
														 // value in structure Parameters.
	std::string step_lengths_f_n;	// f_n is file_name
  Parameters step_numbers;   // Number of steps for every numerical value in
														 //	structure Parameters.
	std::string step_numbers_f_n;	// f_n is file_name
	Parameters orders;         // For each numerical value the is an order in
														 // tree structure of parameters
	std::string orders_f_n; // f_n is file_name


  void Get(std::string & initial_config_file_name,
					 std::string & step_lengths_file_name,
					 std::string & step_numbers_file_name,
					 std::string & orders_file_name);

	void CreateFS(std::string & sim_directory_name);
	void CreateBashForIndependenSim(std::string & bash_file_name,
																	std::string & sim_derectory_name,
																	std::string & prog_file_name,
																	unsigned int proc_num = 1,
																	unsigned int limit_of_time = 240 /*limit time for each simulation (minutes)*/,
																	unsigned int operational_memory_per_proc = 1024,
																	unsigned int memory_per_proc = 0);
	void GetBashForDirInAllCreation(std::string & bash_file_name,
																	std:: string & sim_derectory_name,
																	std::vector<std::string> & vec_file_names);
private:
	void CheckFilesInfo(); //TODO
	std::vector<multi_type> vec_parameter_tuples; // vector of paths and file info
/////auxiliary function for void CreateFS()/////////////////////////////////////
	void WorkWithSubTree(boost::property_tree::ptree & p_tree,
											 std::string s_key,
											 std::string s_last);
	void CreateCatalogesAndFiles(unsigned int tuple_index, std::string relative_path,
															 Parameters parameters); // recursive function
	void ModifyCurrentParameters(Parameters & parameters, int tuple_index, int current_step);
  std::string GetSystemRelevantSlash();
  std::string GetSystemRelevantInit();

  std::string GetSystemRelevantUpperRoot();
////////////////////////////////////////////////////////////////////////////////
	void SetFilesInCataloges();
	//boost::property_tree::ptree FindPropertyWithIndex(int index);
	void CreateDir(std::string & path_to_dir, std::string & dir);
	void BashOpenLoops(std::string & bash_file_name);
	void BashCloseLoops(std::string & bash_file_name);

};


#endif // FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED
