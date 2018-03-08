#include "conf_structures.h"

//initializes parameters
void ParameterStructClass::getFromFile(std::string & file_name) {
    //read parameters from file with name file_name
    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(file_name, pt);

    basic_sim_parameters.seed = pt.get<int>("basic_sim_parameters.seed");
//    basic_sim_parameters.etalon_result_file_name = pt.get<std::string>
//																("basic_sim_parameters.etalon_result_file_name");
    basic_sim_parameters.etalon_result_catalogue = pt.get<std::string>
										("basic_sim_parameters.etalon_result_catalogue");
    basic_sim_parameters.threshold_simulation_num =
      pt.get<int> ("basic_sim_parameters.threshold_simulation_num");

    basic_sim_parameters.evolution.driver_adv =
      pt.get <float> ("basic_sim_parameters.evolution.driver_adv");
    basic_sim_parameters.evolution.mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.mutation_rate");
    basic_sim_parameters.evolution.driver_mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.driver_mutation_rate");
    basic_sim_parameters.evolution.stop_time_cell_num =
      pt.get <int> ("basic_sim_parameters.evolution.stop_time_cell_num");

    basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name =
      pt.get <std::string> ("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name");
    basic_sim_parameters.probes.balk_cells.needle_diam =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.needle_diam");
    basic_sim_parameters.probes.balk_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.use");

		basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name =
			pt.get <std::string> ("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name");
    basic_sim_parameters.probes.tree_cells.cell_num =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.cell_num");
		basic_sim_parameters.probes.tree_cells.in_balk =
			pt.get <int> ("basic_sim_parameters.probes.tree_cells.in_balk");
    basic_sim_parameters.probes.tree_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.use");
//------------------------------------------------------------------------------
		inference_parameters.algorithm_parameters.final_epsilon =
			pt.get <float> ("inference_parameters.algorithm_parameters.final_epsilon");
		inference_parameters.algorithm_parameters.final_iteration =
			pt.get <int> ("inference_parameters.algorithm_parameters.final_iteration");
		inference_parameters.algorithm_parameters.parameter_population_number =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_number");
		inference_parameters.algorithm_parameters.parameter_population_num_alpha =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_num_alpha");

		inference_parameters.algorithm_parameters.p_acc_min =
			pt.get <float> ("inference_parameters.algorithm_parameters.p_acc_min");

		inference_parameters.basic_settings.current_result_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.current_result_file_name");
		inference_parameters.basic_settings.log_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.log_file_name");


		inference_parameters.comparison_parameters.tree_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.tree_dist_type");
		inference_parameters.comparison_parameters.balk_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.balk_dist_type");
		inference_parameters.comparison_parameters.proportion_of_tree_dist =
		  pt.get <float> ("inference_parameters.comparison_parameters.proportion_of_tree_dist");
		inference_parameters.comparison_parameters.threshold_frac_cell_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_cell_diff");
		inference_parameters.comparison_parameters.threshold_frac_mut_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_mut_diff");
		inference_parameters.comparison_parameters.piece_num =
			pt.get <int> ("inference_parameters.comparison_parameters.piece_num");

		inference_parameters.log_scale_search_space.driver_advantage =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_advantage");
		inference_parameters.log_scale_search_space.driver_mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_mutation_rate");
		inference_parameters.log_scale_search_space.mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.mutation_rate");

			p_tree = pt;
};

void ParameterStructClass::setToFile(std::string & file_name) {
  boost::property_tree::ptree pt;

  pt.add("basic_sim_parameters.seed", basic_sim_parameters.seed);
//  pt.add("basic_sim_parameters.etalon_result_file_name",
//				basic_sim_parameters.etalon_result_file_name);
  pt.add("basic_sim_parameters.etalon_result_catalogue",
				basic_sim_parameters.etalon_result_catalogue);
	pt.add("basic_sim_parameters.threshold_simulation_num",
					basic_sim_parameters.threshold_simulation_num);

  pt.add("basic_sim_parameters.evolution.driver_adv",
				basic_sim_parameters.evolution.driver_adv);
  pt.add("basic_sim_parameters.evolution.mutation_rate",
				basic_sim_parameters.evolution.mutation_rate);
  pt.add("basic_sim_parameters.evolution.driver_mutation_rate",
				basic_sim_parameters.evolution.driver_mutation_rate);
  pt.add("basic_sim_parameters.evolution.stop_time_cell_num",
				basic_sim_parameters.evolution.stop_time_cell_num);


	pt.add("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name",
					basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name);
	pt.add("basic_sim_parameters.probes.balk_cells.needle_diam",
					basic_sim_parameters.probes.balk_cells.needle_diam);
  pt.add("basic_sim_parameters.probes.balk_cells.use",
				  basic_sim_parameters.probes.balk_cells.use);
	pt.add("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name",
					basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name);

	pt.add("basic_sim_parameters.probes.tree_cells.cell_num",
					basic_sim_parameters.probes.tree_cells.cell_num);
  pt.add("basic_sim_parameters.probes.tree_cells.in_balk",
					basic_sim_parameters.probes.tree_cells.in_balk);
  pt.add("basic_sim_parameters.probes.tree_cells.use",
				  basic_sim_parameters.probes.tree_cells.use);
//------------------------------------------------------------------------------

  pt.add("inference_parameters.algorithm_parameters.final_epsilon",
				inference_parameters.algorithm_parameters.final_epsilon);
  pt.add("inference_parameters.algorithm_parameters.final_iteration",
				inference_parameters.algorithm_parameters.final_iteration);
  pt.add("inference_parameters.algorithm_parameters.parameter_population_number",
				inference_parameters.algorithm_parameters.parameter_population_number);
  pt.add("inference_parameters.algorithm_parameters.parameter_population_num_alpha",
				inference_parameters.algorithm_parameters.parameter_population_num_alpha);

  pt.add("inference_parameters.algorithm_parameters.p_acc_min",
				inference_parameters.algorithm_parameters.p_acc_min);

  pt.add("inference_parameters.basic_settings.current_result_file_name",
				inference_parameters.basic_settings.current_result_file_name);
  pt.add("inference_parameters.basic_settings.log_file_name",
				inference_parameters.basic_settings.log_file_name);

  pt.add("inference_parameters.comparison_parameters.tree_dist_type",
				inference_parameters.comparison_parameters.tree_dist_type);
  pt.add("inference_parameters.comparison_parameters.balk_dist_type",
				inference_parameters.comparison_parameters.balk_dist_type);
  pt.add("inference_parameters.comparison_parameters.proportion_of_tree_dist",
				inference_parameters.comparison_parameters.proportion_of_tree_dist);
  pt.add("inference_parameters.comparison_parameters.threshold_frac_cell_diff",
				inference_parameters.comparison_parameters.threshold_frac_cell_diff);
  pt.add("inference_parameters.comparison_parameters.threshold_frac_mut_diff",
				inference_parameters.comparison_parameters.threshold_frac_mut_diff);
  pt.add("inference_parameters.comparison_parameters.piece_num",
				inference_parameters.comparison_parameters.piece_num);


  pt.add("inference_parameters.log_scale_search_space.driver_advantage",
				inference_parameters.log_scale_search_space.driver_advantage);
  pt.add("inference_parameters.log_scale_search_space.driver_mutation_rate",
				inference_parameters.log_scale_search_space.driver_mutation_rate);
  pt.add("inference_parameters.log_scale_search_space.mutation_rate",
				inference_parameters.log_scale_search_space.mutation_rate);

  boost::property_tree::info_parser::write_info(file_name, pt);

};

boost::property_tree::ptree ParameterStructClass::getPTree(){
	boost::property_tree::ptree pt;

  pt.add("basic_sim_parameters.seed", basic_sim_parameters.seed);
//  pt.add("basic_sim_parameters.etalon_result_file_name",
//				basic_sim_parameters.etalon_result_file_name);
  pt.add("basic_sim_parameters.etalon_result_catalogue",
				basic_sim_parameters.etalon_result_catalogue);
	pt.add("basic_sim_parameters.threshold_simulation_num",
					basic_sim_parameters.threshold_simulation_num);

  pt.add("basic_sim_parameters.evolution.driver_adv",
				basic_sim_parameters.evolution.driver_adv);
  pt.add("basic_sim_parameters.evolution.mutation_rate",
				basic_sim_parameters.evolution.mutation_rate);
  pt.add("basic_sim_parameters.evolution.driver_mutation_rate",
				basic_sim_parameters.evolution.driver_mutation_rate);
  pt.add("basic_sim_parameters.evolution.stop_time_cell_num",
				basic_sim_parameters.evolution.stop_time_cell_num);

	pt.add("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name",
					basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name);
	pt.add("basic_sim_parameters.probes.balk_cells.needle_diam",
					basic_sim_parameters.probes.balk_cells.needle_diam);
	pt.add("basic_sim_parameters.probes.balk_cells.use",
					basic_sim_parameters.probes.balk_cells.use);

	pt.add("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name",
					basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name);
	pt.add("basic_sim_parameters.probes.tree_cells.cell_num",
					basic_sim_parameters.probes.tree_cells.cell_num);
	pt.add("basic_sim_parameters.probes.tree_cells.in_balk",
					basic_sim_parameters.probes.tree_cells.in_balk);
	pt.add("basic_sim_parameters.probes.tree_cells.use",
					basic_sim_parameters.probes.tree_cells.use);
//------------------------------------------------------------------------------

  pt.add("inference_parameters.algorithm_parameters.final_epsilon",
				inference_parameters.algorithm_parameters.final_epsilon);
  pt.add("inference_parameters.algorithm_parameters.final_iteration",
				inference_parameters.algorithm_parameters.final_iteration);
  pt.add("inference_parameters.algorithm_parameters.parameter_population_number",
				inference_parameters.algorithm_parameters.parameter_population_number);
  pt.add("inference_parameters.algorithm_parameters.parameter_population_num_alpha",
				inference_parameters.algorithm_parameters.parameter_population_num_alpha);


  pt.add("inference_parameters.algorithm_parameters.p_acc_min",
				inference_parameters.algorithm_parameters.p_acc_min);

  pt.add("inference_parameters.basic_settings.current_result_file_name",
				inference_parameters.basic_settings.current_result_file_name);
  pt.add("inference_parameters.basic_settings.log_file_name",
				inference_parameters.basic_settings.log_file_name);


  pt.add("inference_parameters.comparison_parameters.tree_dist_type",
				inference_parameters.comparison_parameters.tree_dist_type);
  pt.add("inference_parameters.comparison_parameters.balk_dist_type",
				inference_parameters.comparison_parameters.balk_dist_type);
  pt.add("inference_parameters.comparison_parameters.proportion_of_tree_dist",
				inference_parameters.comparison_parameters.proportion_of_tree_dist);
  pt.add("inference_parameters.comparison_parameters.threshold_frac_cell_diff",
				inference_parameters.comparison_parameters.threshold_frac_cell_diff);
  pt.add("inference_parameters.comparison_parameters.threshold_frac_mut_diff",
				inference_parameters.comparison_parameters.threshold_frac_mut_diff);
  pt.add("inference_parameters.comparison_parameters.piece_num",
				inference_parameters.comparison_parameters.piece_num);

  pt.add("inference_parameters.log_scale_search_space.driver_advantage",
				inference_parameters.log_scale_search_space.driver_advantage);
  pt.add("inference_parameters.log_scale_search_space.driver_mutation_rate",
				inference_parameters.log_scale_search_space.driver_mutation_rate);
  pt.add("inference_parameters.log_scale_search_space.mutation_rate",
				inference_parameters.log_scale_search_space.mutation_rate);

	return pt;
};

void ParameterStructClass::getFromPTree(boost::property_tree::ptree pt) {

    basic_sim_parameters.seed = pt.get<int>("basic_sim_parameters.seed");
//    basic_sim_parameters.etalon_result_file_name = pt.get<std::string>
//																("basic_sim_parameters.etalon_result_file_name");
    basic_sim_parameters.etalon_result_catalogue = pt.get<std::string>
										("basic_sim_parameters.etalon_result_catalogue");
    basic_sim_parameters.threshold_simulation_num =
      pt.get<int> ("basic_sim_parameters.threshold_simulation_num");

    basic_sim_parameters.evolution.driver_adv =
      pt.get <float> ("basic_sim_parameters.evolution.driver_adv");
    basic_sim_parameters.evolution.mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.mutation_rate");
    basic_sim_parameters.evolution.driver_mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.driver_mutation_rate");
    basic_sim_parameters.evolution.stop_time_cell_num =
      pt.get <int> ("basic_sim_parameters.evolution.stop_time_cell_num");

    basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name =
      pt.get <std::string> ("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name");
    basic_sim_parameters.probes.balk_cells.needle_diam =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.needle_diam");
    basic_sim_parameters.probes.balk_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.use");

		basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name =
			pt.get <std::string> ("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name");
    basic_sim_parameters.probes.tree_cells.cell_num =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.cell_num");
		basic_sim_parameters.probes.tree_cells.in_balk =
			pt.get <int> ("basic_sim_parameters.probes.tree_cells.in_balk");
    basic_sim_parameters.probes.tree_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.use");
//------------------------------------------------------------------------------
		inference_parameters.algorithm_parameters.final_epsilon =
			pt.get <float> ("inference_parameters.algorithm_parameters.final_epsilon");
		inference_parameters.algorithm_parameters.final_iteration =
			pt.get <int> ("inference_parameters.algorithm_parameters.final_iteration");
		inference_parameters.algorithm_parameters.parameter_population_number =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_number");
		inference_parameters.algorithm_parameters.parameter_population_num_alpha =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_num_alpha");


		inference_parameters.algorithm_parameters.p_acc_min =
			pt.get <float> ("inference_parameters.algorithm_parameters.p_acc_min");

		inference_parameters.basic_settings.current_result_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.current_result_file_name");
		inference_parameters.basic_settings.log_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.log_file_name");


		inference_parameters.comparison_parameters.tree_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.tree_dist_type");
		inference_parameters.comparison_parameters.balk_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.balk_dist_type");
		inference_parameters.comparison_parameters.proportion_of_tree_dist =
		  pt.get <float> ("inference_parameters.comparison_parameters.proportion_of_tree_dist");
		inference_parameters.comparison_parameters.threshold_frac_cell_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_cell_diff");
		inference_parameters.comparison_parameters.threshold_frac_mut_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_mut_diff");
		inference_parameters.comparison_parameters.piece_num =
			pt.get <int> ("inference_parameters.comparison_parameters.piece_num");


		inference_parameters.log_scale_search_space.driver_advantage =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_advantage");
		inference_parameters.log_scale_search_space.driver_mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_mutation_rate");
		inference_parameters.log_scale_search_space.mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.mutation_rate");

			p_tree = pt;

};

void ParameterStructClass::printToConsole() {

    std::cout<< "basic_sim_parameters.seed  = " <<
			basic_sim_parameters.seed <<std::endl;
//    std::cout<< "basic_sim_parameters.etalon_result_file_name  = " <<
//			basic_sim_parameters.etalon_result_file_name <<std::endl;
    std::cout<< "basic_sim_parameters.etalon_result_catalogue = " <<
			basic_sim_parameters.etalon_result_catalogue <<std::endl;
    std::cout<< "basic_sim_parameters.threshold_simulation_num  = " <<
			basic_sim_parameters.threshold_simulation_num <<std::endl;

    std::cout<< "basic_sim_parameters.evolution.driver_adv  = " <<
			basic_sim_parameters.evolution.driver_adv <<std::endl;
		std::cout<< "basic_sim_parameters.evolution.mutation_rate  = " <<
			basic_sim_parameters.evolution.mutation_rate <<std::endl;
    std::cout<< "basic_sim_parameters.evolution.driver_mutation_rate = " <<
			basic_sim_parameters.evolution.driver_mutation_rate <<std::endl;
    std::cout<< "basic_sim_parameters.evolution.stop_time_cell_num = " <<
			basic_sim_parameters.evolution.stop_time_cell_num <<std::endl;



	std::cout<<"basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name" <<
					basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name<<std::endl;
	std::cout<<"basic_sim_parameters.probes.balk_cells.needle_diam" <<
					basic_sim_parameters.probes.balk_cells.needle_diam<<std::endl;
	std::cout<<"basic_sim_parameters.probes.balk_cells.use" <<
					basic_sim_parameters.probes.balk_cells.use<<std::endl;

	std::cout<<"basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name" <<
					basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name<<std::endl;
	std::cout<<"basic_sim_parameters.probes.tree_cells.cell_num" <<
					basic_sim_parameters.probes.tree_cells.cell_num<<std::endl;
	std::cout<<"basic_sim_parameters.probes.tree_cells.in_balk" <<
					basic_sim_parameters.probes.tree_cells.in_balk<<std::endl;
	std::cout<<"basic_sim_parameters.probes.tree_cells.use" <<
					basic_sim_parameters.probes.tree_cells.use<<std::endl;
//------------------------------------------------------------------------------
		std::cout<< "inference_parameters.algorithm_parameters.final_epsilon = " <<
			inference_parameters.algorithm_parameters.final_epsilon <<std::endl;
		std::cout<< "inference_parameters.algorithm_parameters.final_iteration = " <<
			inference_parameters.algorithm_parameters.final_iteration <<std::endl;
		std::cout<< "inference_parameters.algorithm_parameters.parameter_population_number = " <<
			inference_parameters.algorithm_parameters.parameter_population_number <<std::endl;
		std::cout<< "inference_parameters.algorithm_parameters.parameter_population_num_alpha = " <<
			inference_parameters.algorithm_parameters.parameter_population_num_alpha <<std::endl;


		std::cout<< "inference_parameters.algorithm_parameters.p_acc_min = " <<
			inference_parameters.algorithm_parameters.p_acc_min <<std::endl;
		std::cout<< "inference_parameters.basic_settings.current_result_file_name = " <<
			inference_parameters.basic_settings.current_result_file_name <<std::endl;
		std::cout<< "inference_parameters.basic_settings.log_file_name = " <<
			inference_parameters.basic_settings.log_file_name <<std::endl;


		std::cout<< "inference_parameters.comparison_parameters.tree_dist_type = " <<
			inference_parameters.comparison_parameters.tree_dist_type <<std::endl;
		std::cout<< "inference_parameters.comparison_parameters.balk_dist_type = " <<
			inference_parameters.comparison_parameters.balk_dist_type <<std::endl;
		std::cout<< "inference_parameters.comparison_parameters.proportion_of_tree_dist" <<
				inference_parameters.comparison_parameters.proportion_of_tree_dist <<std::endl;

		std::cout<< "inference_parameters.comparison_parameters.threshold_frac_cell_diff = " <<
			inference_parameters.comparison_parameters.threshold_frac_cell_diff <<std::endl;
		std::cout<< "inference_parameters.comparison_parameters.threshold_frac_mut_diff = " <<
			inference_parameters.comparison_parameters.threshold_frac_mut_diff <<std::endl;
		std::cout<< "inference_parameters.comparison_parameters.piece_num = "	<<
			inference_parameters.comparison_parameters.piece_num << std::endl;


		std::cout<< "inference_parameters.log_scale_search_space.driver_advantage = " <<
			inference_parameters.log_scale_search_space.driver_advantage <<std::endl;
		std::cout<< "inference_parameters.log_scale_search_space.driver_mutation_rate = " <<
			inference_parameters.log_scale_search_space.driver_mutation_rate <<std::endl;
		std::cout<< "inference_parameters.log_scale_search_space.mutation_rate = " <<
			inference_parameters.log_scale_search_space.mutation_rate <<std::endl;
};

std::string ParameterStructClass::getSearchSpaceParsString() {

    std::string search_space_pars_string;
    search_space_pars_string =
				boost::lexical_cast<std::string>(basic_sim_parameters.evolution.driver_adv) + " "
			+ boost::lexical_cast<std::string>(basic_sim_parameters.evolution.mutation_rate) + " "
			+ boost::lexical_cast<std::string>(basic_sim_parameters.evolution.driver_mutation_rate);
    return search_space_pars_string;
};
////////////////////////////////////////////////////////////////////////////////

void getEtalonPartParameters(BasicSimParameters & basic_sim_parameters){
    boost::property_tree::ptree pt;
    std::string full_par_file_name = basic_sim_parameters.related_path_to_par_file +
																		 basic_sim_parameters.par_file_name;
    boost::property_tree::info_parser::read_info(
					full_par_file_name, pt);

    basic_sim_parameters.seed = pt.get<int> ("basic_sim_parameters.seed");
		basic_sim_parameters.etalon_result_catalogue = pt.get<std::string>
										("basic_sim_parameters.etalon_result_catalogue");
//    basic_sim_parameters.etalon_result_file_name =
//			pt.get<std::string> ("basic_sim_parameters.etalon_result_file_name");

    basic_sim_parameters.threshold_simulation_num =
      pt.get<int> ("basic_sim_parameters.threshold_simulation_num");
    basic_sim_parameters.evolution.driver_adv =
      pt.get <float> ("basic_sim_parameters.evolution.driver_adv");
    basic_sim_parameters.evolution.mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.mutation_rate");
    basic_sim_parameters.evolution.driver_mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.driver_mutation_rate");
    basic_sim_parameters.evolution.stop_time_cell_num =
      pt.get <int> ("basic_sim_parameters.evolution.stop_time_cell_num");

		basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name =
      pt.get <std::string> ("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name");
    basic_sim_parameters.probes.balk_cells.needle_diam =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.needle_diam");
    basic_sim_parameters.probes.balk_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.use");

		basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name =
			pt.get <std::string> ("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name");
    basic_sim_parameters.probes.tree_cells.cell_num =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.cell_num");
		basic_sim_parameters.probes.tree_cells.in_balk =
			pt.get <int> ("basic_sim_parameters.probes.tree_cells.in_balk");
    basic_sim_parameters.probes.tree_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.use");

};

void getInferencePartParameters(BasicSimParameters & basic_sim_parameters,
																InferenceParameters & inference_parameters) {
    boost::property_tree::ptree pt;
    std::string full_par_file_name = basic_sim_parameters.related_path_to_par_file +
																		 basic_sim_parameters.par_file_name;
    boost::property_tree::info_parser::read_info(
					full_par_file_name, pt);
//We take only part of parameters due to external info

//    basic_sim_parameters.etalon_result_file_name =
//			pt.get<std::string> ("basic_sim_parameters.etalon_result_file_name");
		basic_sim_parameters.etalon_result_catalogue =
			pt.get<std::string> ("basic_sim_parameters.etalon_result_catalogue");
    basic_sim_parameters.threshold_simulation_num =
      pt.get<int> ("basic_sim_parameters.threshold_simulation_num");

    basic_sim_parameters.evolution.stop_time_cell_num =
      pt.get <int> ("basic_sim_parameters.evolution.stop_time_cell_num");

    basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name =
      pt.get <std::string> ("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name");
    basic_sim_parameters.probes.balk_cells.needle_diam =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.needle_diam");
    basic_sim_parameters.probes.balk_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.use");

		basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name =
			pt.get <std::string> ("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name");
    basic_sim_parameters.probes.tree_cells.cell_num =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.cell_num");
		basic_sim_parameters.probes.tree_cells.in_balk =
			pt.get <int> ("basic_sim_parameters.probes.tree_cells.in_balk");
    basic_sim_parameters.probes.tree_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.use");

		inference_parameters.comparison_parameters.piece_num =
			pt.get <int> ("inference_parameters.comparison_parameters.piece_num");
		inference_parameters.comparison_parameters.tree_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.tree_dist_type");
		inference_parameters.comparison_parameters.balk_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.balk_dist_type");
		inference_parameters.comparison_parameters.proportion_of_tree_dist =
		  pt.get <float> ("inference_parameters.comparison_parameters.proportion_of_tree_dist");

		inference_parameters.comparison_parameters.threshold_frac_cell_diff =
		  pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_cell_diff");
		inference_parameters.comparison_parameters.threshold_frac_mut_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_mut_diff");
};

void getAllParameters(BasicSimParameters & basic_sim_parameters,
											InferenceParameters & inference_parameters) {
		std::string full_par_file_name = basic_sim_parameters.related_path_to_par_file +
																		 basic_sim_parameters.par_file_name;
    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(
					full_par_file_name, pt);
    basic_sim_parameters.seed = pt.get<int> ("basic_sim_parameters.seed");
//    basic_sim_parameters.etalon_result_file_name =
//			pt.get<std::string> ("basic_sim_parameters.etalon_result_file_name");
		basic_sim_parameters.etalon_result_catalogue =
			pt.get<std::string>("basic_sim_parameters.etalon_result_catalogue");
    basic_sim_parameters.threshold_simulation_num =
      pt.get<int> ("basic_sim_parameters.threshold_simulation_num");
    basic_sim_parameters.evolution.driver_adv =
      pt.get <float> ("basic_sim_parameters.evolution.driver_adv");
    basic_sim_parameters.evolution.mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.mutation_rate");
    basic_sim_parameters.evolution.driver_mutation_rate =
      pt.get <float>("basic_sim_parameters.evolution.driver_mutation_rate");
    basic_sim_parameters.evolution.stop_time_cell_num =
      pt.get <int> ("basic_sim_parameters.evolution.stop_time_cell_num");

    basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name =
      pt.get <std::string> ("basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name");
    basic_sim_parameters.probes.balk_cells.needle_diam =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.needle_diam");
    basic_sim_parameters.probes.balk_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.balk_cells.use");

		basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name =
			pt.get <std::string> ("basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name");
    basic_sim_parameters.probes.tree_cells.cell_num =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.cell_num");
		basic_sim_parameters.probes.tree_cells.in_balk =
			pt.get <int> ("basic_sim_parameters.probes.tree_cells.in_balk");
    basic_sim_parameters.probes.tree_cells.use =
      pt.get <int> ("basic_sim_parameters.probes.tree_cells.use");
//------------------------------------------------------------------------------
		inference_parameters.algorithm_parameters.final_epsilon =
			pt.get <float> ("inference_parameters.algorithm_parameters.final_epsilon");
		inference_parameters.algorithm_parameters.final_iteration =
			pt.get <int> ("inference_parameters.algorithm_parameters.final_iteration");
		inference_parameters.algorithm_parameters.parameter_population_number =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_number");
		inference_parameters.algorithm_parameters.parameter_population_num_alpha =
			pt.get <int> ("inference_parameters.algorithm_parameters.parameter_population_num_alpha");
		inference_parameters.algorithm_parameters.p_acc_min =
			pt.get <float> ("inference_parameters.algorithm_parameters.p_acc_min");

		inference_parameters.basic_settings.current_result_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.current_result_file_name");
		inference_parameters.basic_settings.log_file_name =
			pt.get <std::string> ("inference_parameters.basic_settings.log_file_name");

		inference_parameters.comparison_parameters.tree_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.tree_dist_type");
		inference_parameters.comparison_parameters.balk_dist_type =
			pt.get <int> ("inference_parameters.comparison_parameters.balk_dist_type");
		inference_parameters.comparison_parameters.proportion_of_tree_dist =
		  pt.get <float> ("inference_parameters.comparison_parameters.proportion_of_tree_dist");
		inference_parameters.comparison_parameters.threshold_frac_cell_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_cell_diff");
		inference_parameters.comparison_parameters.threshold_frac_mut_diff =
			pt.get <float> ("inference_parameters.comparison_parameters.threshold_frac_mut_diff");
		inference_parameters.comparison_parameters.piece_num =
			pt.get <int> ("inference_parameters.comparison_parameters.piece_num");

		inference_parameters.log_scale_search_space.driver_advantage =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_advantage");
		inference_parameters.log_scale_search_space.driver_mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.driver_mutation_rate");
		inference_parameters.log_scale_search_space.mutation_rate =
			pt.get <float> ("inference_parameters.log_scale_search_space.mutation_rate");
};


std::string getSystemRelevantInit(){
#if defined __linux
	std::string slash("./");
#elif defined __APPLE__
  std::string slash("./");
#else
	std::string slash("");
#endif
	return slash;
};

std::string getSystemRelevantOneUpperRoot(){
	std::string slash = getSystemRelevantSlash();
	std::string two_points("..") ;
	std::string upper_catalogue(two_points + slash);
	return upper_catalogue;
};

std::string getSystemRelevantTwoUpperRoot() {
	std::string slash = getSystemRelevantSlash();
	std::string two_points("..") ;
	std::string upper_catalogue(two_points + slash + two_points + slash);
	return upper_catalogue;
};

std::string getSystemRelevantSlash(){
#if defined __linux
	std::string slash("/");
#elif defined __APPLE__
  std::string slash("/");
#else
	std::string slash("\\");
#endif
	return slash;
};

std::string getProgName(){
	std::string prog_name("cancer_3d.exe");
	return prog_name;
};

std::string getABCProgName(){
	std::string prog_name("multi.exe");
	return prog_name;
};

std::string getCatalogueProgName(){
	std::string catalogue("sim_prog");
	return catalogue;
};

std::string getCatalogueABCProgName(){
	std::string catalogue("multi");
	return catalogue;
};



void err(char *reason)
{
  std::cout <<reason<<std::endl ;
#ifdef __WIN32
  system("pause");
#endif
  exit(0) ;
};

void err(const char *reason)
{
  std::cout <<reason<<std::endl ;
#ifdef __WIN32
  system("pause") ;
#endif
  exit(0) ;
}

bool fileExists(const char *name) {
    std::ifstream my_file(name);
    if(my_file.fail()){ //File does not exist
        return false;
    } else { //otherwise, file exists
        return true;
    };
};

bool fileRemoveIfExists(const char *name){
		if(fileExists(name)){
					remove(name);
					return true;
		};
		return false;
};
