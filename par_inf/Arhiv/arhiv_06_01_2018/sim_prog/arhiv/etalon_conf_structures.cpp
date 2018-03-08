#include "etalon_conf_structures.h"

//initializes parameters
void Parameters::GetFromFile(std::string & file_name_) {
    //read parameters from file with name file_name
    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(file_name_, pt);

    seed = pt.get<int> ("parameters.seed");
    file_name = pt.get<std::string> ("parameters.file_name");
    catalog_file_name = pt.get<std::string> ("parameters.catalog_file_name");
    threshold_simulation_num =
      pt.get<int> ("parameters.threshold_simulation_num");
    evolution.driver_adv =
      pt.get <float> ("parameters.evolution.driver_adv");
    evolution.mutation_rate =
      pt.get <float>("parameters.evolution.mutation_rate");
    evolution.driver_mutation_rate =
      pt.get <float>("parameters.evolution.driver_mutation_rate");
    evolution.stop_time_cell_num =
      pt.get <int> ("parameters.evolution.stop_time_cell_num");

    probe.cell_num =
      pt.get <int> ("parameters.probe.cell_num");
    probe.is_random =
			pt.get <int> ("parameters.probe.is_random");

    probe.delta_area.del_x =
      pt.get <int> ("parameters.probe.delta_area.del_x");
    probe.delta_area.del_y =
      pt.get <int> ("parameters.probe.delta_area.del_y");
    probe.delta_area.del_z =
      pt.get <int> ("parameters.probe.delta_area.del_z");
			p_tree = pt;
};

void Parameters::SetToFile(std::string & file_name_) {
  boost::property_tree::ptree pt;

  pt.add("parameters.seed", seed);
  pt.add("parameters.file_name", file_name);
  pt.add("parameters.catalog_file_name", catalog_file_name);
	pt.add("parameters.threshold_simulation_num", threshold_simulation_num);

  pt.add("parameters.evolution.driver_adv", evolution.driver_adv);
  pt.add("parameters.evolution.mutation_rate", evolution.mutation_rate);
  pt.add("parameters.evolution.driver_mutation_rate", evolution.driver_mutation_rate);
  pt.add("parameters.evolution.stop_time_cell_num", evolution.stop_time_cell_num);

  pt.add("parameters.probe.is_random", probe.is_random);
  pt.add("parameters.probe.cell_num", probe.cell_num);
  pt.add("parameters.probe.delta_area.del_x", probe.delta_area.del_x);
  pt.add("parameters.probe.delta_area.del_y", probe.delta_area.del_y);
  pt.add("parameters.probe.delta_area.del_z", probe.delta_area.del_z);

  boost::property_tree::info_parser::write_info(file_name_, pt);


};

boost::property_tree::ptree Parameters::GetPTree() {
	boost::property_tree::ptree pt;

  pt.add("parameters.seed", seed);
  pt.add("parameters.file_name", file_name);
  pt.add("parameters.catalog_file_name", catalog_file_name);
	pt.add("parameters.threshold_simulation_num", threshold_simulation_num);

  pt.add("parameters.evolution.driver_adv", evolution.driver_adv);
  pt.add("parameters.evolution.mutation_rate", evolution.mutation_rate);
  pt.add("parameters.evolution.driver_mutation_rate", evolution.driver_mutation_rate);
  pt.add("parameters.evolution.stop_time_cell_num", evolution.stop_time_cell_num);

  pt.add("parameters.probe.cell_num", probe.cell_num);
	pt.add("parameters.probe.is_random", probe.is_random);
  pt.add("parameters.probe.delta_area.del_x", probe.delta_area.del_x);
  pt.add("parameters.probe.delta_area.del_y", probe.delta_area.del_y);
  pt.add("parameters.probe.delta_area.del_z", probe.delta_area.del_z);

	return pt;
};

void Parameters::GetFromPTree(boost::property_tree::ptree pt){

    seed = pt.get<int> ("parameters.seed");
    file_name = pt.get<std::string> ("parameters.file_name");
    catalog_file_name = pt.get<std::string> ("parameters.catalog_file_name");
    threshold_simulation_num =
      pt.get<int> ("parameters.threshold_simulation_num");

    evolution.driver_adv =
      pt.get <float> ("parameters.evolution.driver_adv");
    evolution.mutation_rate =
      pt.get <float>("parameters.evolution.mutation_rate");
    evolution.driver_mutation_rate =
      pt.get <float>("parameters.evolution.driver_mutation_rate");
    evolution.stop_time_cell_num =
      pt.get <int> ("parameters.evolution.stop_time_cell_num");

    probe.cell_num =
      pt.get <int> ("parameters.probe.cell_num");
    probe.is_random =
			pt.get <int> ("parameters.probe.is_random");

    probe.delta_area.del_x =
      pt.get <int> ("parameters.probe.delta_area.del_x");
    probe.delta_area.del_y =
      pt.get <int> ("parameters.probe.delta_area.del_y");
    probe.delta_area.del_z =
      pt.get <int> ("parameters.probe.delta_area.del_z");

			p_tree = pt;

};

void Parameters::PrintToConsole() {

		std::cout<< "parameters.seed = " <<
								seed <<std::endl;
		std::cout<< "parameters.file_name = " <<
								file_name <<std::endl;
		std::cout<< "parameters.catalog_file_name = " <<
								 catalog_file_name <<std::endl;
		std::cout<< "parameters.threshold_simulation_num = " <<
								 threshold_simulation_num <<std::endl;

		std::cout<< "parameters.evolution.driver_adv = " <<
								 evolution.driver_adv<< std::endl;
		std::cout<< "parameters.evolution.mutation_rate = " <<
								 evolution.mutation_rate << std::endl;
		std::cout<< "parameters.evolution.driver_mutation_rate = " <<
								 evolution.driver_mutation_rate << std::endl;
		std::cout<< "parameters.evolution.stop_time_cell_num = " <<
								 evolution.stop_time_cell_num<< std::endl;

		std::cout<< "parameters.probe.cell_num = " <<
								 probe.cell_num << std::endl;
		std::cout<< "parameters.probe.is_random = " <<
							   probe.is_random << std::endl;

		std::cout<< "parameters.probe.delta_area.del_x = " <<
								 probe.delta_area.del_x << std::endl;
		std::cout<< "parameters.probe.delta_area.del_y = " <<
								 probe.delta_area.del_y << std::endl;
		std::cout<< "parameters.probe.delta_area.del_z = " <<
								 probe.delta_area.del_z << std::endl;
};

