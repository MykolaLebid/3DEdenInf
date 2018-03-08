#include "prog_parameters.h"
//
//
std::string GetFullProgName(){
	std::string full_prog_name(GetProgName());
	full_prog_name = GetSystemRelevantOneUpperRoot() + GetCatalogueProgName()
	+ GetSystemRelevantSlash() + full_prog_name;
	return full_prog_name;
};

int GetRandNumber(Vars & vars){
		std::random_device rd;
		std::uniform_int_distribution<int> rnd(1,2000000000);
    return rnd(vars.gen);
};

void SimulateAndWrite2File(const ParPoint par_point,
													 const AllSimParameters & conf_settings, const int rand) {
    std::string command_line;
    command_line = GetSystemRelevantInit() + GetFullProgName() + " "
								 + conf_settings.basic_sim_parameters.related_path_to_par_file + " "
								 + conf_settings.basic_sim_parameters.par_file_name + " "
								 + boost::lexical_cast<std::string>(par_point.driver_adv) + " "
								 + boost::lexical_cast<std::string>(par_point.mut_rate) + " "
								 + boost::lexical_cast<std::string>(par_point.driver_mut_rate) + " "
								 + boost::lexical_cast<std::string>(rand);
		//std::cout << command_line << std::endl;
    //int i;
    //std::cin >> i ;
    system (command_line.c_str());
};

int GetNumLinesInComparisonFile(const AllSimParameters & conf_settings) {
	std::string full_comparison_file_name = GetFullComparisonFileName(conf_settings);
	std::ifstream file_tree_comparison;
	file_tree_comparison.open(full_comparison_file_name.c_str());

	int number_of_lines = 0;
	std::string line;
	while (std::getline(file_tree_comparison, line)) number_of_lines++;
	file_tree_comparison.close();

	return number_of_lines;
};

////void GetTreeComparisonVecFromFile(const Parameters & pars,
////													        std::vector<TreeComparison> & tree_comparison_vec) {
////	std::string full_comparison_file_name = GetFullComparisonFileName(pars);
////	std::ifstream file_tree_comparison;
////
////	file_tree_comparison.open(full_comparison_file_name.c_str());
////	if (file_tree_comparison.is_open()==false) {
////		std::cout<<"err..."<<std::endl;
////	} else {
////		par_pop_num = pars.conf_settings.inference_parameters.
////											algorithm_parameters.parameter_population_number;
////		for(int i = 0; i < par_pop_num; i++) {
////			TreeComparison A;
////			file_tree_comparison >> A.par_point.driver_adv
////													 >> A.par_point.mut_rate
////													 >> A.par_point.driver_mut_rate
////													 >> A.results.tree_dist
////													 >> A.results.delta_mut
////													 >> A.results.num_mut;
////			if (!file_tree_comparison.eof()) tree_comparison_vec.push_back(A);
////		};
//// };
////}
//
void GetTreeComparisonVecFromFile(const AllSimParameters & conf_settings, const int diff,
														std::vector<TreeComparison> & tree_comparison_vec) {
	std::string full_comparison_file_name = GetFullComparisonFileName(conf_settings);
	std::ifstream file_tree_comparison;

	file_tree_comparison.open(full_comparison_file_name.c_str());

	if (file_tree_comparison.is_open()==false) {
		std::cout<<"err..."<<std::endl;
	} else {
		for(int i = 0; i < diff; i++) {
			TreeComparison A;
			file_tree_comparison >> A.par_point.driver_adv
													 >> A.par_point.mut_rate
													 >> A.par_point.driver_mut_rate
													 >> A.results.balk_dist
													 >> A.results.tree_dist
													 >> A.results.dist;
			if (!file_tree_comparison.eof()) tree_comparison_vec.push_back(A);
		};
 };
}

void SetTreeComparisonVec2File(const std::vector<TreeComparison> & tree_comparison_vec,
															 const AllSimParameters & conf_settings, const int number_of_iteration) {
	std::ofstream file_tree_comparison;
	std::string file_name;
	file_name = GetPathToEtalonResultCatalogue(conf_settings) + GetSystemRelevantSlash() +
						boost::lexical_cast<std::string>(number_of_iteration) + "_" +
						conf_settings.inference_parameters.basic_settings.current_result_file_name;
	file_tree_comparison.open(file_name.c_str());

	if (file_tree_comparison.is_open()==false) std::cout<<"err..."<<std::endl;

	file_tree_comparison.precision(5);
	file_tree_comparison.setf(std::ios::fixed, std::ios::floatfield);

	for(unsigned int i = 0; i < tree_comparison_vec.size(); i++)
		file_tree_comparison <<
			tree_comparison_vec[i].par_point.driver_adv << "\t" <<
			tree_comparison_vec[i].par_point.mut_rate << "\t" <<
			tree_comparison_vec[i].par_point.driver_mut_rate << "\t" <<
			tree_comparison_vec[i].results.tree_dist << "\t" <<
			tree_comparison_vec[i].results.weight << std::endl;
	    file_tree_comparison.close();
};

//// Paths and full file names
//////////////////////////////////////////////////////////////////////////////////
std::string GetPathToEtalonResultCatalogue(const AllSimParameters & conf_settings) {
	std::string file_name =
			conf_settings.basic_sim_parameters.related_path_to_par_file +
			conf_settings.basic_sim_parameters.etalon_result_catalogue +
			GetSystemRelevantSlash();
	return file_name;
};

std::string GetFullComparisonFileName(const AllSimParameters & conf_settings) {
	std::string file_name = GetPathToEtalonResultCatalogue(conf_settings);
	file_name += conf_settings.inference_parameters.basic_settings.
								current_result_file_name;
	return file_name;
};

void AddToLogFile(const AllSimParameters & conf_settings,
									const unsigned int iteration_num,
									const unsigned int trials,
									const TreeComparison & relative_error,
									const TreeComparison & sample_var){
	std::ofstream log_tree_dist_comparison;

	std::string full_log_file_name = GetPathToEtalonResultCatalogue(conf_settings) +
					conf_settings.inference_parameters.basic_settings.log_file_name;

	log_tree_dist_comparison.open(full_log_file_name.c_str(), std::ios_base::app);
	if (log_tree_dist_comparison.is_open() == false) std::cout<<"err..."<<std::endl;
	log_tree_dist_comparison.precision(5);
	log_tree_dist_comparison.setf(std::ios::fixed, std:: ios::floatfield);
	log_tree_dist_comparison << iteration_num                             <<"\t"<<
														        trials                     	        <<"\t"<<
		conf_settings.basic_sim_parameters.evolution.driver_adv             <<"\t"<<
	                            relative_error.par_point.driver_adv       <<"\t"<<
	                            sample_var.par_point.driver_adv           <<"\t"<<

	  conf_settings.basic_sim_parameters.evolution.mutation_rate          <<"\t"<<
	                            relative_error.par_point.mut_rate         <<"\t"<<
	                            sample_var.par_point.mut_rate             <<"\t"<<
	  conf_settings.basic_sim_parameters.evolution.driver_mutation_rate   <<"\t"<<
	                            relative_error.par_point.driver_mut_rate  <<"\t"<<
	                            sample_var.par_point.driver_mut_rate<< std::endl;
	log_tree_dist_comparison.close();
};

//////////////////////////////////////////////////////////////////////////////////
//
//////initializes parameters
////void GetPars(Parameters & pars) {
////  boost::property_tree::ptree pt;
////  boost::property_tree::info_parser::read_info(pars.par_file_name, pt);
////  // begin* Name a directory with etalon_mut_vec.dat (etalon mutation vector)
////  pars.catalog_file_name = pt.get<string>
////		("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file");
////  pars.comparison_file_name = pt.get<string>
////    ("Parameters.ProbeComparisonParameters.comparison_file_name");
////  pars.number_of_parameters = 3;
////	pars.algo.final_epsilon = pt.get<float>
////	 ("algorithm_parameters.final_epsilon");
////	pars.algo.final_iteration = pt.get<float>
////	 ("algorithm_parameters.final_iteration");
////  pars.vars.max_dist = 1000;
////  pars.algo.parameter_population_num = pt.get<float>
////		("algorithm_parameters.parameter_population_number");
////  pars.algo.parameter_population_num_alpha =
////		pars.algo.parameter_population_num/2;
////  pars.algo.p_acc_min = pt.get<float>
////	 ("algorithm_parameters.p_acc_min");
////
////  pars.search_space.driver_adv.min_ = pt.get<float>
////		("search_space.driver_advantage.min");
////  pars.search_space.driver_adv.max_ = pt.get<float>
////		("search_space.driver_advantage.max");
////  pars.search_space.driver_adv.true_val = pt.get<float>
////		("search_space.driver_advantage.true");
////
////  pars.search_space.mut_rate.min_ = pt.get<float>
////		("search_space.mutation_rate.min");
////  pars.search_space.mut_rate.max_ =  pt.get<float>
////		("search_space.mutation_rate.max");
////  pars.search_space.mut_rate.true_val = pt.get<float>
////		("search_space.mutation_rate.true");
////
////  pars.search_space.driver_mut_rate.min_ = pt.get<float>
////		("search_space.driver_mutation_rate.min");
////  pars.search_space.driver_mut_rate.max_ = pt.get<float>
////		("search_space.driver_mutation_rate.max");
////  pars.search_space.driver_mut_rate.true_val = pt.get<float>
////    ("search_space.driver_mutation_rate.true");
////
////	float stoping_level = 0.0001;
////  float driver_adv_range = pars.search_space.driver_adv.max_-
////                          pars.search_space.driver_adv.min_;
////  float mut_rate_range = pars.search_space.mut_rate.max_-
////                          pars.search_space.mut_rate.min_;
////  float driver_mut_rate_range = pars.search_space.driver_mut_rate.max_-
////                          pars.search_space.driver_mut_rate.min_;
////
////  if (driver_adv_range/pars.search_space.driver_adv.min_ <= stoping_level) {
////		pars.search_space.driver_adv.is_mutable = false;
////	} else {
////		pars.search_space.driver_adv.is_mutable = true;
////	}
////  if (mut_rate_range/pars.search_space.mut_rate.min_ <= stoping_level) {
////		pars.search_space.mut_rate.is_mutable = false;
////  } else {
////    pars.search_space.mut_rate.is_mutable = true;
////  };
////  if (driver_mut_rate_range/pars.search_space.driver_mut_rate.min_ <= stoping_level) {
////		pars.search_space.driver_mut_rate.is_mutable = false;
////  } else {
////  	pars.search_space.driver_mut_rate.is_mutable =true;
////  };
////
////  pars.vars.p_acc = 1000;
////  pars.log_file_name = "logs.dat";
////};


float GetLeftBorderDriverAdvantage(const AllSimParameters & conf_settings){
		int log_driver_advantage = conf_settings.inference_parameters.
															   log_scale_search_space.driver_advantage;
		return (conf_settings.basic_sim_parameters.evolution.driver_adv /
						  std::pow(10.0, log_driver_advantage));
};

float GetRigthBorderDriverAdvantage(const AllSimParameters & conf_settings){
		int log_driver_advantage = conf_settings.inference_parameters.
																 log_scale_search_space.driver_advantage;
		return (conf_settings.basic_sim_parameters.evolution.driver_adv *
						  std::pow(10.0,log_driver_advantage));
}

float GetLeftBorderMutationRate(const AllSimParameters & conf_settings){
	int log_mutation_rate = conf_settings.inference_parameters.
														log_scale_search_space.mutation_rate;
	return (conf_settings.basic_sim_parameters.evolution.mutation_rate /
					  std::pow(10.0,log_mutation_rate));
};

float GetRigthBorderMutationRate(const AllSimParameters & conf_settings){
	int log_mutation_rate = conf_settings.inference_parameters.
														log_scale_search_space.mutation_rate;
	return (conf_settings.basic_sim_parameters.evolution.mutation_rate *
				    std::pow(10.0,log_mutation_rate));
};

float GetLeftBorderDriverMutationRate(const AllSimParameters & conf_settings){
	int log_driver_mut_rate = conf_settings.inference_parameters.
															log_scale_search_space.driver_mutation_rate;

	return  (conf_settings.basic_sim_parameters.evolution.driver_mutation_rate /
						 std::pow(10.0,log_driver_mut_rate));
};

float GetRigthBorderDriverMutationRate(const AllSimParameters & conf_settings){
	int log_driver_mut_rate =
			conf_settings.inference_parameters.log_scale_search_space.driver_mutation_rate;
	return (conf_settings.basic_sim_parameters.evolution.driver_mutation_rate *
					  std::pow(10.0,log_driver_mut_rate));
};

bool CompareTwoTrees(const TreeComparison tree_results_1, const TreeComparison tree_results_2) {
    return (tree_results_1.results.tree_dist <
					  tree_results_2.results.tree_dist);
};

int GetSearchParNumber(){return 3;};
