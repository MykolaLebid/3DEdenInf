#ifndef FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED
#define FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED

#include <conf_structures.h>
//#include <boost/lexical_cast.hpp>

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
  ParameterStructClass initial_config; // Structure with all initial configurations.
														 //	This structure will be as a template for all
														 // config files.
	std::string initial_config_f_n; // f_n is file_name
  ParameterStructClass step_lengths;   // Step length for every numerical
														 // value in ParameterStructClass.
	std::string step_lengths_f_n;	// f_n is file_name
  ParameterStructClass step_numbers;   // Number of steps for every numerical value in
														 //	structure ParameterStructClass.
	std::string step_numbers_f_n;	// f_n is file_name
	ParameterStructClass orders; // For each numerical value the is an order in
														 // tree structure of ParameterStructClass
	std::string orders_f_n; // f_n is file_name

	std::string upper_dir_with_config_files;



  void Get(std::string & initial_config_file_name,
					 std::string & step_lengths_file_name,
					 std::string & step_numbers_file_name,
					 std::string & orders_file_name);

	void CreateFS(std::string & sim_directory_name);
	void CreateBashForIndependenSim(const std::string & bash_file_name,
																	const std::string & sim_directory_name,
																	const std::string & initial_config_f_n,
																	unsigned int proc_num = 1,
																	unsigned int limit_of_time = 240 /*limit time for each simulation (minutes)*/,
																	unsigned int operational_memory_per_proc = 1024,
																	unsigned int memory_per_proc = 0);

	void CreateBashForABCInference(const std::string & bash_file_name_ABC,
																 const std::string & sim_directory_name,
																 const std::string & initial_config_f_n,
																 unsigned int proc_num = 10,
																 unsigned int limit_of_time = 240 /*limit time for each simulation (minutes)*/,
																 unsigned int operational_memory_per_proc = 1024,
																 unsigned int memory_per_proc = 0);

	void GetBashForDirInAllCreation(std::string & bash_file_name,
																	std::string & sim_directory_name,
																	std::vector<std::string> & vec_file_names);
	std::vector<std::string> GetPathVector(std::string init_dir);
private:
	void CheckFilesInfo(); //TODO
	std::vector<multi_type> vec_parameter_tuples; // vector of paths and file info
/////auxiliary function for void CreateFS()/////////////////////////////////////
	void WorkWithSubTree(boost::property_tree::ptree & p_tree,
											 std::string s_key,
											 std::string s_last);
	void CreateCataloguesAndFiles(unsigned int tuple_index, std::string relative_path,
															  ParameterStructClass parameters); // recursive function
	void ModifyCurrentParameters(ParameterStructClass & parameters,
															 int tuple_index, int current_step);
//  std::string GetSystemRelevantSlash();
//  std::string GetSystemRelevantInit();

  std::string GetSystemRelevantUpperRoot();
////////////////////////////////////////////////////////////////////////////////
	void SetFilesInCatalogues();
	//boost::property_tree::ptree FindPropertyWithIndex(int index);
	void CreateDir(std::string & path_to_dir, std::string & dir);
	void WriteIntroduction();
	void BashOpenLoops(const std::string & bash_file_name);
	void BashCloseLoops(const std::string & bash_file_name);
	void BashHead(const std::string & bash_file_name, const unsigned int proc_num);

	void GetCatalogueNames(unsigned int tuple_index, std::string relative_path,
												 std::vector<std::string> & path_vec);
	// supplementary function for GetPathVector();
};


#endif // FILE_SYS_FOR_PARALLEL_COMPUTATIONS_H_INCLUDED
