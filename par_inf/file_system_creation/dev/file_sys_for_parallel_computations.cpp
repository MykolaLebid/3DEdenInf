#include "file_sys_for_parallel_computations.h"



//#include<iostream>
////---------------conf file---------------------------begin
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/info_parser.hpp>
////---------------conf file-----------------------------end
//
////initializes parameters
void StatisticsFS::get(std::string & initial_config_file_name,
											 std::string & step_lengths_file_name,
											 std::string & step_numbers_file_name,
											 std::string & orders_file_name) {
  std::string dir_path = ".." + getSystemRelevantSlash();

  std::string local_name = dir_path + initial_config_file_name;
	initial_config.getFromFile(local_name);
	initial_config_f_n = initial_config_file_name;

	local_name = dir_path + step_lengths_file_name;
	step_lengths.getFromFile(local_name);
	step_lengths_f_n = step_lengths_file_name;

	local_name = dir_path + step_numbers_file_name;
	step_numbers.getFromFile(local_name);
	step_numbers_f_n = step_numbers_file_name;

	local_name = dir_path + orders_file_name;
	orders.getFromFile(local_name);
	orders_f_n = orders_file_name;

	checkFilesInfo(); //TODO
	//boost::filesystem::path path_to_file(".\\Development");
	//std::cout <<   boost::filesystem::file_size(path_to_file) << '\n';
};
//creates directories and setting files according to 4 main files
void StatisticsFS::createFS(std::string & sim_directory_name){
	//take all information from files recursively and write to vec_parameter_tuples

	boost::property_tree::ptree pt_initial_config = initial_config.getPTree();
  std::string st = "";
	workWithSubTree(pt_initial_config, st, st);
	// sort  vec_parameter_tuples by orders
	std::sort(vec_parameter_tuples.begin(), vec_parameter_tuples.end(),
						[](const multi_type& a, const multi_type& b)
            {
                return std::get<4>(a) < std::get<4>(b);
						} );
	std::string full_sim_directory_name	=	getSystemRelevantTwoUpperRoot() +
																				sim_directory_name;

//	const boost::filesystem::path path;
//	boost::filesystem::create_directory( path );
//
//	if ( boost::filesystem::exists( path ) )
//	{
////        boost::filesystem::remove_all( path );
//	};
///////////////////////Create upper directory///////////////////////////////////
  std::string command_line ("mkdir ");
	command_line += full_sim_directory_name;
	system(command_line.c_str());
////////////////////////////////////////////////////////////////////////////////
	createCataloguesAndFiles(0, sim_directory_name, initial_config);
//  for(auto x: vec_parameter_tuples) {
//		std::cout<< std::get<0>(x) << "  " <<
//								std::get<1>(x) << "  " <<
//								std::get<2>(x) << "  " <<
//								std::get<3>(x) << "  " <<
//								std::get<4>(x) << "  " <<std::endl;
//  }
};

void StatisticsFS::initFSforStatistics(std::string & sim_directory_name){
	//take all information from files recursively and write to vec_parameter_tuples

	boost::property_tree::ptree pt_initial_config = initial_config.getPTree();
  std::string st = "";
	workWithSubTree(pt_initial_config, st, st);
	// sort  vec_parameter_tuples by orders
	std::sort(vec_parameter_tuples.begin(), vec_parameter_tuples.end(),
						[](const multi_type& a, const multi_type& b)
            {
                return std::get<4>(a) < std::get<4>(b);
						} );
	std::string full_sim_directory_name	=	getSystemRelevantTwoUpperRoot() +
																				sim_directory_name;
};

//std::string StatisticsFS::GetSystemRelevantTwoUpperRoot(){
//	std::string slash = GetSystemRelevantSlash();
//	std::string two_points("..") ;
//	std::string upper_catalogue(two_points + slash + two_points + slash);
//	return upper_catalogue;
//};

//std::string StatisticsFS::GetSystemRelevantSlash(){
//
//#if defined __linux
//	std::string slash("/");
//#elif defined __APPLE__
//  std::string slash("/");
//#else
//	std::string slash("\\");
//#endif
//	return slash;
//};
//



void StatisticsFS::createCataloguesAndFiles(unsigned int tuple_index,
																						std::string relative_path,
																						ParameterStructClass parameters){
	if (tuple_index < vec_parameter_tuples.size()) {
		std::string slash = getSystemRelevantSlash();
		std::string underscore("_") ;

		std::string catalogue_name = relative_path + slash;
		int step_number = std::get<3>(vec_parameter_tuples.at(tuple_index));
		std::string command_line ("mkdir ");

		for (int i = 0; i < step_number; i++) {

			const std::string local_catalogue_name = catalogue_name +
				                std::get<1>(vec_parameter_tuples.at(tuple_index)) +
				                underscore +
				                boost::lexical_cast<std::string>(i);
			std::string current_command_line = command_line +
																				 getSystemRelevantTwoUpperRoot() +
																				 local_catalogue_name;
			system(current_command_line.c_str());

			modifyCurrentParameters(parameters, tuple_index, i);
			if ((tuple_index + 1) == vec_parameter_tuples.size()) {
				std::string path_to_conf_file = getSystemRelevantTwoUpperRoot() + local_catalogue_name +
																				slash + initial_config_f_n;
				parameters.setToFile(path_to_conf_file);
			};
			createCataloguesAndFiles(tuple_index + 1, local_catalogue_name, parameters);
		};
	};
};

void StatisticsFS::modifyCurrentParameters(ParameterStructClass & parameters,
																					 int tuple_index,
																					 int current_step){

			boost::property_tree::ptree p_tree = parameters.getPTree();
			boost::property_tree::ptree p_initial_config_tree = initial_config.getPTree();
			boost::property_tree::ptree p_step_lengths_tree = step_lengths.getPTree();

			std::string path_to_par = std::get<0>(vec_parameter_tuples.at(tuple_index));
			float init_var = p_initial_config_tree.get<float>(path_to_par, 0);
		  float step_length = p_step_lengths_tree.get<float>(path_to_par, 0);
		  float new_var = init_var + step_length * current_step;
			p_tree.put(path_to_par, new_var);

 		  parameters.getFromPTree(p_tree);
};

void StatisticsFS::workWithSubTree(boost::property_tree::ptree & p_tree,
																	 std::string s_key,
																	 std::string s_last){
	if (p_tree.empty()) {
//		std::cout<<s_key << std::endl;
//    std::cout<<s_last << std::endl;
    boost::property_tree::ptree pt_orders = orders.getPTree();
    //orders.PrintToConsole();

    if ( pt_orders.get <int> (s_key) != 0) {
				 boost::property_tree::ptree pt_step_lengths = step_lengths.getPTree();
				 boost::property_tree::ptree pt_step_numbers = step_numbers.getPTree();
				 boost::property_tree::ptree pt_orders = orders.getPTree();
				 multi_type local_tuple =
				 std::make_tuple(s_key,
												 s_last,
												 pt_step_lengths.get <float> (s_key),
												 pt_step_numbers.get <int> (s_key),
												 pt_orders.get <int> (s_key));
				 vec_parameter_tuples.push_back(local_tuple);
    };

  } else {
		for (boost::property_tree::ptree::iterator pos = p_tree.begin(); pos != p_tree.end();) {
      std::string new_s_key;
      if (s_key!="") {
				new_s_key = s_key +"."+ pos->first;
      } else {
      	new_s_key = pos->first;
			};
      workWithSubTree(pos->second, new_s_key, pos->first);
      ++pos;
		};

  };
};


void StatisticsFS::bashOpenLoops(const std::string & bash_file_name){
		std::ofstream bash_file;
	  bash_file.open(bash_file_name, std::ios_base::app);
		unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");

		for(unsigned int i = 0; i < number_of_parameters; i++) {
			bash_file <<locale_space<<"for((i_"<<i<<"=0;i_"<<i<<
								"<="<<(std::get<3>(vec_parameter_tuples.at(i))-1)<<
								";i_"<<i<<"++))"<<std::endl;
			bash_file <<locale_space<< "do" << std::endl;
			locale_space += space;
		};
		bash_file.close();
};

//void BashOpenLoopsForSearchSpace(const std::string & bash_file_name,
//																 const int number_of_parameters_in_search_space,
//																 const int number_repeats){
//		std::ofstream bash_file;
//	  bash_file.open(bash_file_name, std::ios_base::app);
//		std::string space("  ");
//		std::string locale_space("");
//		for(unsigned int j = 0; j < number_of_parameters_in_search_space; j++) {
//			bash_file <<locale_space<<"for((j_"<<j<<"=0;j_"<<j<<
//								"<"<<number<<
//								";j_"<<j<<"++))"<<std::endl;
//			bash_file <<locale_space<< "do" << std::endl;
//			locale_space += space;
//		};
//		bash_file.close();
//
//
//};

void StatisticsFS::bashCloseLoops(const std::string & bash_file_name){
		std::ofstream bash_file;
	  bash_file.open(bash_file_name, std::ios_base::app);
		unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");
		for(unsigned int i = 0; i < number_of_parameters; i++) {
			locale_space = "";
			for(unsigned int j = number_of_parameters - (i+1) ; j > 0; j--)
				locale_space += space;
			bash_file << locale_space << "done" << std::endl;
			locale_space += space;
		};
		bash_file.close();
};

void StatisticsFS::bashHead(const std::string & bash_file_name,
														const unsigned int proc_num){
		std::ofstream bash_file;
	  bash_file.open(bash_file_name, std::ios_base::app);
		bash_file <<  "# set number of processors for each thread" << std::endl <<
									"export OMP_NUM_THREADS = " << proc_num<< std::endl;

		bash_file <<
		"# send all jobs to Euler" << std::endl <<
		"# -n = number of cores (should be the same as OMP_NUM_THREADS)" << std::endl <<
		"# -W = limit time for each thread"<< std::endl <<
		"# -R rusage[mem=XXX] - operational memory per processor" <<  std::endl <<
		"# -R rusage[scratch=XXX] - disk memory per processor" <<  std::endl;
		bash_file.close();
};

//void StatisticsFS::CreateBashForTheSameSimInAllCataloges();
void StatisticsFS::createBashEtalonSims( const std::string & bash_file_name,
																				 const std::string & sim_directory_name,
																				 const std::string & initial_config_f_n,
																				 unsigned int proc_num,
																				 unsigned int limit_of_time,
																				 unsigned int operational_memory_per_proc,
																				 unsigned int memory_per_proc){
		fileRemoveIfExists(bash_file_name.c_str());
		bashHead(bash_file_name, proc_num);
		bashOpenLoops(bash_file_name);
    unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");

		for(unsigned int j = number_of_parameters; j > 0; j--)
			locale_space += space;

		std::ofstream bash_file;

	  bash_file.open(bash_file_name, std::ios_base::app);
		bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
							<< limit_of_time << " "
							<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
							<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
							<< getSystemRelevantInit() /*<< GetSystemRelevantOneUpperRoot()
							<< GetCatalogueABCProgName() << GetSystemRelevantSlash()*/ << getProgName() <<  " "
							<< getLocaleNamePathToSettingFile(sim_directory_name) << " "
							<< initial_config_f_n <<std::endl;

		bashCloseLoops(bash_file_name);
		bash_file.close();


};

//void StatisticsFS::CreateBashForProbeSimViaMainSettingFile(
//																	const std::string & bash_file_name,
//																	const std::string & sim_directory_name,
//																	const std::string & initial_config_f_n,
//																	unsigned int proc_num,
//																	unsigned int limit_of_time/*limit time for each simulation (minutes)*/,
//																	unsigned int operational_memory_per_proc,
//																	unsigned int memory_per_proc){
//		FileRemoveIfExists(bash_file_name.c_str());
//		BashHead(bash_file_name, proc_num);
//		BashOpenLoops(bash_file_name);
//    unsigned int number_of_parameters = vec_parameter_tuples.size();
//		std::string space("  ");
//		std::string locale_space("");
//
//		for(unsigned int j = number_of_parameters; j > 0; j--)
//			locale_space += space;
//
//		std::ofstream bash_file;
//
//	  bash_file.open(bash_file_name, std::ios_base::app);
//		bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
//							<< limit_of_time << " "
//							<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
//							<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
//							<< GetSystemRelevantInit() << GetProgName() <<  " "
//							<< GetSystemRelevantOneUpperRoot() << " "
//							<< initial_config_f_n <<std::endl;
//
//		BashCloseLoops(bash_file_name);
//		bash_file.close();
//
//
//
//};

void StatisticsFS::createBashCompareDist(   const std::string & bash_file_name_compare_dist,
																						const std::string & sim_directory_name,
																						const std::string & init_conf_file_name,
																						const unsigned int point_number,
																						const unsigned int iteration_number,
																						const float left_scale,
																						const float right_scale,
																						const unsigned int iteration_mode,
																						unsigned int proc_num,
																						unsigned int limit_of_time,
															              unsigned int operational_memory_per_proc,
															              unsigned int memory_per_proc){
	fileRemoveIfExists(bash_file_name_compare_dist.c_str());
	bashHead(bash_file_name_compare_dist, proc_num);
	bashOpenLoops(bash_file_name_compare_dist);

	unsigned int number_of_parameters = vec_parameter_tuples.size();
	std::string space("  ");
	std::string locale_space("");

	for(unsigned int j = number_of_parameters; j > 0; j--)
		locale_space += space;

	std::ofstream bash_file;
	bash_file.open(bash_file_name_compare_dist, std::ios_base::app);
	bash_file.precision(5);
	bash_file.setf(std::ios::fixed, std:: ios::floatfield);
	bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
						<< limit_of_time << " "
						<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
						<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
					  << getSystemRelevantInit() /*<< GetSystemRelevantOneUpperRoot()
						<< GetCatalogueProgName() << GetSystemRelevantSlash()*/ << getProgName() <<  " "
						<< getLocaleNamePathToSettingFile(sim_directory_name) << " "
						<< initial_config_f_n << " " ;

	for (unsigned int i = 0; i < number_of_parameters; i++)
		bash_file <<"$i_"<<boost::lexical_cast<std::string>(i);

	bash_file << " " << point_number << " " << iteration_number<< " "
	          << left_scale<< " " << right_scale << " "  << iteration_mode <<  std::endl;
	bashCloseLoops(bash_file_name_compare_dist);
	bash_file.close();
};

void StatisticsFS::createBashStatistics( const std::string & bash_file_name_statistics,
																				 const std::string & sim_directory_name,
																				 const std::string & init_conf_file_name,
																				 const std::string & statistics_type){
	fileRemoveIfExists(bash_file_name_statistics.c_str());
	//BashHead(bash_file_name_statistics, proc_num);
	bashOpenLoops(bash_file_name_statistics);

	unsigned int number_of_parameters = vec_parameter_tuples.size();
	std::string space("  ");
	std::string locale_space("");

	for(unsigned int j = number_of_parameters; j > 0; j--)
		locale_space += space;

	std::ofstream bash_file;
	bash_file.open(bash_file_name_statistics, std::ios_base::app);
	bash_file.precision(5);
	bash_file.setf(std::ios::fixed, std:: ios::floatfield);
	bash_file << locale_space << getSystemRelevantInit() << getProgName() <<  " "
						<< getLocaleNamePathToSettingFile(sim_directory_name) << " "
						<< initial_config_f_n << " " ;

  bash_file << statistics_type <<  std::endl;
	bashCloseLoops(bash_file_name_statistics);
	bash_file.close();
};

std::string StatisticsFS::getLocaleNamePathToSettingFile(const std::string & sim_directory_name){
	std::string locale_name_of_setting_file("");
	locale_name_of_setting_file += getSystemRelevantOneUpperRoot() +
																 sim_directory_name +
																 getSystemRelevantSlash();
  unsigned int number_of_parameters = vec_parameter_tuples.size();

	for (unsigned int i = 0; i < number_of_parameters; i++)
	locale_name_of_setting_file += std::get<1>(vec_parameter_tuples.at(i)) +
																 "_$i_" + boost::lexical_cast<std::string>(i) +
																 getSystemRelevantSlash();
	return locale_name_of_setting_file;
};

void StatisticsFS::createBashForABCInference(const std::string & bash_file_name_ABC,
															 const std::string & sim_directory_name,
															 const std::string & initial_config_f_n,
															 unsigned int proc_num,
															 unsigned int limit_of_time,
															 unsigned int operational_memory_per_proc,
															 unsigned int memory_per_proc){
		fileRemoveIfExists(bash_file_name_ABC.c_str());
		bashHead(bash_file_name_ABC, proc_num);
		bashOpenLoops(bash_file_name_ABC);
    unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");

		for(unsigned int j = number_of_parameters; j > 0; j--)
			locale_space += space;




		std::ofstream bash_file;
	  bash_file.open(bash_file_name_ABC, std::ios_base::app);
		bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
							<< limit_of_time << " "
							<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
							<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
							<< getSystemRelevantInit() /*<< GetSystemRelevantOneUpperRoot()
							<< GetCatalogueProgName() << GetSystemRelevantSlash()*/ << getABCProgName() <<  " "
							<< getLocaleNamePathToSettingFile(sim_directory_name) << " "
							<< initial_config_f_n << " ";

		for (unsigned int i = 0; i < number_of_parameters; i++)
			bash_file <<"$i_"<<boost::lexical_cast<std::string>(i);
		bash_file << std::endl;

		bashCloseLoops(bash_file_name_ABC);
		bash_file.close();

};



void StatisticsFS::getBashForDirInAllCreation(std::string & bash_file_name,
																							std::string & sim_directory_name,
																              std::vector<std::string> & vec_of_dirs) {
	fileRemoveIfExists(bash_file_name.c_str());
	bashOpenLoops(bash_file_name);
	std::ofstream bash_file;
	bash_file.open(bash_file_name, std::ios_base::app);
  unsigned int number_of_parameters = vec_parameter_tuples.size();
	std::string space("  ");
	std::string locale_space("");
	for(unsigned int j = number_of_parameters; j > 0; j--)
			locale_space += space;

	for_each (vec_of_dirs.begin(), vec_of_dirs.end(), [&](std::string str) {
								std::string locale_name_of_setting_file("");
								for (unsigned int i = 0; i < number_of_parameters; i++)
								locale_name_of_setting_file += std::get<1>(vec_parameter_tuples.at(i)) +
																				"_$i_"+ boost::lexical_cast<std::string>(i) +
																				getSystemRelevantSlash();
								bash_file << locale_space<< "mkdir -p "
													<< sim_directory_name << getSystemRelevantSlash()
													<< locale_name_of_setting_file << str <<std::endl;
						});

	bashCloseLoops(bash_file_name);

	bash_file.close();
};


std::vector<std::string> StatisticsFS::getPathVector(std::string init_dir) {
	std::vector<std::string> path_vec;
	getCatalogueNames(0, init_dir, path_vec);
	return path_vec;
};

void StatisticsFS::getCatalogueNames(unsigned int tuple_index,
																		 std::string relative_path,
																		 std::vector<std::string> & path_vec) {
	unsigned int vec_size = vec_parameter_tuples.size();
	if (tuple_index < vec_size) {
		std::string slash = getSystemRelevantSlash();
		std::string underscore("_") ;

		std::string catalogue_name = relative_path + slash;
		int step_number = std::get<3>(vec_parameter_tuples.at(tuple_index));
		for (int i = 0; i < step_number; i++) {

			const std::string local_catalogue_name = catalogue_name +
				                std::get<1>(vec_parameter_tuples.at(tuple_index)) +
				                underscore +
				                boost::lexical_cast<std::string>(i);

			if ((tuple_index + 1) == vec_size) {
				std::string path_to_conf_file = getSystemRelevantTwoUpperRoot() +
																				local_catalogue_name;
				path_vec.push_back(path_to_conf_file);
			} else {
				getCatalogueNames(tuple_index + 1,
						local_catalogue_name, path_vec);
			};
		};
	};
};
//TODO
void StatisticsFS::checkFilesInfo(){

};


