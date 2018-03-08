#include "file_sys_for_parallel_computations.h"



//#include<iostream>
////---------------conf file---------------------------begin
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/info_parser.hpp>
////---------------conf file-----------------------------end
//
////initializes parameters
void StatisticsFS::Get(std::string & initial_config_file_name,
											 std::string & step_lengths_file_name,
											 std::string & step_numbers_file_name,
											 std::string & orders_file_name) {
  std::string dir_path = ".." + GetSystemRelevantSlash();

  std::string local_name = dir_path + initial_config_file_name;
	initial_config.GetFromFile(local_name);
	initial_config_f_n = initial_config_file_name;

	local_name = dir_path + step_lengths_file_name;
	step_lengths.GetFromFile(local_name);
	step_lengths_f_n = step_lengths_file_name;

	local_name = dir_path + step_numbers_file_name;
	step_numbers.GetFromFile(local_name);
	step_numbers_f_n = step_numbers_file_name;

	local_name = dir_path + orders_file_name;
	orders.GetFromFile(local_name);
	orders_f_n = orders_file_name;

	CheckFilesInfo(); //TODO
	//boost::filesystem::path path_to_file(".\\Development");
	//std::cout <<   boost::filesystem::file_size(path_to_file) << '\n';
};

void StatisticsFS::CreateFS(std::string & sim_directory_name){
	//take all information from files recursively and write to vec_parameter_tuples

	boost::property_tree::ptree pt_initial_config = initial_config.getPTree();
  std::string st = "";
	WorkWithSubTree(pt_initial_config, st, st);
	// sort  vec_parameter_tuples by orders
	std::sort(vec_parameter_tuples.begin(), vec_parameter_tuples.end(),
						[](const multi_type& a, const multi_type& b)
            {
                return std::get<4>(a) < std::get<4>(b);
						} );
	std::string full_sim_directory_name	=	GetSystemRelevantTwoUpperRoot() +
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
	CreateCataloguesAndFiles(0, sim_directory_name, initial_config);
//  for(auto x: vec_parameter_tuples) {
//		std::cout<< std::get<0>(x) << "  " <<
//								std::get<1>(x) << "  " <<
//								std::get<2>(x) << "  " <<
//								std::get<3>(x) << "  " <<
//								std::get<4>(x) << "  " <<std::endl;
//  }
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


void StatisticsFS::CreateCataloguesAndFiles(unsigned int tuple_index,
																					 std::string relative_path,
																					 ParameterStructClass parameters){
	if (tuple_index < vec_parameter_tuples.size()) {
		std::string slash = GetSystemRelevantSlash();
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
																				 GetSystemRelevantTwoUpperRoot() +
																				 local_catalogue_name;
			system(current_command_line.c_str());
			ModifyCurrentParameters(parameters, tuple_index, i);
			if ((tuple_index + 1) == vec_parameter_tuples.size()) {
				std::string path_to_conf_file = GetSystemRelevantTwoUpperRoot() + local_catalogue_name +
																				slash + initial_config_f_n;

				parameters.SetToFile(path_to_conf_file);
			};
			CreateCataloguesAndFiles(tuple_index + 1, local_catalogue_name, parameters);
		};
	};
};

void StatisticsFS::ModifyCurrentParameters(ParameterStructClass & parameters,
																					 int tuple_index,
																					 int current_step){

			boost::property_tree::ptree p_tree = parameters.getPTree();
			boost::property_tree::ptree p_initial_config_tree = initial_config.getPTree();
			boost::property_tree::ptree p_step_lengths_tree = step_lengths.getPTree();
			//boost::property_tree::ptree p_step_numbers_tree = step_numbers.GetPTree();

			std::string path_to_par = std::get<0>(vec_parameter_tuples.at(tuple_index));
			float init_var = p_initial_config_tree.get<float>(path_to_par, 0);
		  float step_length = p_step_lengths_tree.get<float>(path_to_par, 0);
		  float new_var = init_var + step_length * current_step;
			p_tree.put(path_to_par, new_var);

 		  parameters.GetFromPTree(p_tree);
};

void StatisticsFS::WorkWithSubTree(boost::property_tree::ptree & p_tree,
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
      WorkWithSubTree(pos->second, new_s_key, pos->first);
      ++pos;
		};

  };
};


void StatisticsFS::BashOpenLoops(const std::string & bash_file_name) {
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

void StatisticsFS::BashCloseLoops(const std::string & bash_file_name){
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

void StatisticsFS::BashHead(const std::string & bash_file_name,
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

void StatisticsFS::CreateBashForIndependenSim(const std::string & bash_file_name,
																							const std::string & sim_derectory_name,
																	            const std::string & initial_config_f_n,
																	            unsigned int proc_num,
																	            unsigned int limit_of_time,
																	            unsigned int operational_memory_per_proc,
																	            unsigned int memory_per_proc) {
		BashHead(bash_file_name, proc_num);
		BashOpenLoops(bash_file_name);
    unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");
		std::string locale_name_of_setting_file("");
		for(unsigned int j = number_of_parameters; j > 0; j--)
			locale_space += space;

    locale_name_of_setting_file += GetSystemRelevantOneUpperRoot() +
																	 sim_derectory_name +
																	 GetSystemRelevantSlash();

		for (unsigned int i = 0; i < number_of_parameters; i++)
			locale_name_of_setting_file += std::get<1>(vec_parameter_tuples.at(i)) +
																		 "_$i_" + boost::lexical_cast<std::string>(i) +
																		 GetSystemRelevantSlash();
		std::ofstream bash_file;
	  bash_file.open(bash_file_name, std::ios_base::app);
		bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
							<< limit_of_time << " "
							<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
							<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
							<< GetSystemRelevantInit() /*<< GetSystemRelevantOneUpperRoot()
							<< GetCatalogueABCProgName() << GetSystemRelevantSlash()*/ << GetProgName() <<  " "
							<< locale_name_of_setting_file << " "
							<< initial_config_f_n <<std::endl;

		BashCloseLoops(bash_file_name);
		bash_file.close();


};

void StatisticsFS::CreateBashForABCInference(const std::string & bash_file_name_ABC,
															 const std::string & sim_directory_name,
															 const std::string & initial_config_f_n,
															 unsigned int proc_num,
															 unsigned int limit_of_time,
															 unsigned int operational_memory_per_proc,
															 unsigned int memory_per_proc){
		BashHead(bash_file_name_ABC, proc_num);
		BashOpenLoops(bash_file_name_ABC);
    unsigned int number_of_parameters = vec_parameter_tuples.size();
		std::string space("  ");
		std::string locale_space("");
		std::string locale_name_of_setting_file("");
		for(unsigned int j = number_of_parameters; j > 0; j--)
			locale_space += space;

    locale_name_of_setting_file += GetSystemRelevantOneUpperRoot() +
																	 sim_directory_name +
																	 GetSystemRelevantSlash();

		for (unsigned int i = 0; i < number_of_parameters; i++)
		locale_name_of_setting_file += std::get<1>(vec_parameter_tuples.at(i)) +
																	   "_$i_" + boost::lexical_cast<std::string>(i) +
																		 GetSystemRelevantSlash();
		std::ofstream bash_file;
	  bash_file.open(bash_file_name_ABC, std::ios_base::app);
		bash_file << locale_space <<"bsub -n "<< proc_num <<" -W "
							<< limit_of_time << " "
							<< "-R \"rusage[mem="<< operational_memory_per_proc <<"]\" "
							<< " "  << "-R \"rusage[scratch="<< memory_per_proc <<"]\" "
							<< GetSystemRelevantInit() /*<< GetSystemRelevantOneUpperRoot()
							<< GetCatalogueProgName() << GetSystemRelevantSlash()*/ << GetABCProgName() <<  " "
							<< locale_name_of_setting_file << " "
							<< initial_config_f_n << " ";

		for (unsigned int i = 0; i < number_of_parameters; i++)
			bash_file <<"$i_"<<boost::lexical_cast<std::string>(i);
		bash_file <<std::endl;

		BashCloseLoops(bash_file_name_ABC);
		bash_file.close();

};



void StatisticsFS::GetBashForDirInAllCreation(std::string & bash_file_name,
																							std::string & sim_directory_name,
																              std::vector<std::string> & vec_of_dirs) {
	BashOpenLoops(bash_file_name);
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
																				GetSystemRelevantSlash();
								bash_file << locale_space<< "mkdir -p "
													<< sim_directory_name << GetSystemRelevantSlash()
													<< locale_name_of_setting_file << str <<std::endl;
						});

	BashCloseLoops(bash_file_name);

	bash_file.close();
};


std::vector<std::string> StatisticsFS::GetPathVector(std::string init_dir) {
	std::vector<std::string> path_vec;
	GetCatalogueNames(0, init_dir, path_vec);
	return path_vec;
};

void StatisticsFS::GetCatalogueNames(unsigned int tuple_index,
																		 std::string relative_path,
																		 std::vector<std::string> & path_vec) {
	unsigned int vec_size = vec_parameter_tuples.size();
	if (tuple_index < vec_size) {
		std::string slash = GetSystemRelevantSlash();
		std::string underscore("_") ;

		std::string catalogue_name = relative_path + slash;
		int step_number = std::get<3>(vec_parameter_tuples.at(tuple_index));
		for (int i = 0; i < step_number; i++) {

			const std::string local_catalogue_name = catalogue_name +
				                std::get<1>(vec_parameter_tuples.at(tuple_index)) +
				                underscore +
				                boost::lexical_cast<std::string>(i);

			if ((tuple_index + 1) == vec_size) {
				std::string path_to_conf_file = GetSystemRelevantTwoUpperRoot() +
																				local_catalogue_name;
				path_vec.push_back(path_to_conf_file);
			} else {
				GetCatalogueNames(tuple_index + 1,
						local_catalogue_name, path_vec);
			};
		};
	};
};
//TODO
void StatisticsFS::CheckFilesInfo(){

};


