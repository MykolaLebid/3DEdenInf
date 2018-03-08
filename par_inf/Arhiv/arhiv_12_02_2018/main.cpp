// Copyright 2017 Lebid Mykola all rights reserved

#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "params.h"
#include "classes.h"
#include "sample_processing.h"

////////////////////////////////////////////////////
//TODO delete
//#define __MAIN
//#include <stdio.h>
//#include <stdlib.h>
//#include <cmath>
//#include <stdexcept>// out_of_range
//#include <fstream>
//checking existence of a file
//#include <boost/filesystem.hpp>
//timer
//#include <chrono>
////////////////////////////////////////////////////

#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL)
  #error no method defined!
#endif

#if (defined(VON_NEUMANN_NEIGHBOURHOOD) && defined(MOORE_NEIGHBOURHOOD))
  #error both VON_NEUMANN_NEIGHBOURHOOD and MOORE_NEIGHBOURHOOD defined!
#endif

#if (!defined(VON_NEUMANN_NEIGHBOURHOOD) && !defined(MOORE_NEIGHBOURHOOD))
  #error neither VON_NEUMANN_NEIGHBOURHOOD nor MOORE_NEIGHBOURHOOD defined!
#endif

//extern char *NUM;
//extern const char *NUM_E;
//extern std::string STR;
//extern int /*RAND, sample, treatment, max_size,*/ nsam;
//extern double tt ;//time

//int max_size;
float driver_adv;
float driver_prob;
float gama;
float death0, growth0 = 1;

int num_adv;

//int sample=0;

//search_min_max
void initMinMaxBordars(BasicSimParameters & pars){
//Mykola search of min and max for all cells in all lesions
	int minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;

  for (auto & cell:cells){
    if (cell.x<minx) minx=cell.x;
    if (cell.x>maxx) maxx=cell.x;
    if (cell.y<miny) miny=cell.y;
    if (cell.y>maxy) maxy=cell.y;
    if (cell.z<minz) minz=cell.z;
    if (cell.z>maxz) maxz=cell.z;
  }

  maxx++;
  maxy++;
	maxz++;

	pars.probes.population_borders.max_x = maxx;
	pars.probes.population_borders.max_y = maxy;
  pars.probes.population_borders.max_z = maxz;
	pars.probes.population_borders.min_x = minx;
	pars.probes.population_borders.min_y = miny;
  pars.probes.population_borders.min_z = minz;
};

bool run(BasicSimParameters & pars, bool is_etalon_sim) {
   std::string path = pars.related_path_to_par_file +
											pars.etalon_result_catalogue;
   init(is_etalon_sim, path);
   int s = 0;
   reset(pars.evolution.driver_adv, pars.evolution.driver_mutation_rate);
   while ((main_proc(pars.evolution.stop_time_cell_num,
                     2,-1, -1)==1) &&
          (s < pars.threshold_simulation_num)) {
        s++ ; reset(pars.evolution.driver_adv, pars.evolution.driver_mutation_rate);
    } // initial growth until max size is reached
	 initMinMaxBordars(pars);
	 if (s < pars.threshold_simulation_num) return true;
	 else return false;

};


void initSimGlobalValues(const BasicSimParameters & pars) {
	//RAND				= pars.seed;
	//max_size 		= pars.evolution.stop_time_cell_num;
	driver_adv  = pars.evolution.driver_adv;
	gama  			= pars.evolution.mutation_rate;
	driver_prob = pars.evolution.driver_mutation_rate;
	growth0 		= 1;
	death0 			= growth0 * (1 - pars.evolution.driver_adv);
	_srand48(pars.seed);

	//Name a directory with etalon_mut_vec.dat (etalon mutation vector)
	//std::string str_1 = pars.etalon_result_file_name;
	//std::string str_1 = pars.etalon_result_file_name;
	//char *cstr = new char[str_1.length() + 1];
	//strcpy(cstr, str_1.c_str());
  //NUM = cstr;
};

void InitParamSearch (char *name) {
    //read parameters from file "model_info.info"
    std::string config = std::string(name);
    //cout<<config<<endl;

    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(config, pt);

    //RAND = pt.get<int> ("basic_settings.seed");
    //Name a directory with etalon_mut_vec.dat (etalon mutation vector)
    //std::string str_1 = pt.get<std::string> ("basic_settings.file_name");
    //char *cstr = new char[str_1.length() + 1];
    //strcpy(cstr, str_1.c_str());
    //NUM = cstr;


    // do stuff
    //cout<< "NUM" << NUM << endl;
    //cout<< "str.c_str()="<<str.c_str()<<endl;
    //cout<< "(char*) str.c_str()="<<(char*) str.c_str()<<endl;
    //max_size    =  pt.get <float>("evolution_settings.stopping_time");
    driver_adv  =  pt.get <float>("search_space.driver_advantage.true");
    gama        =  pt.get <float>("search_space.mutation_rate.true");
    driver_prob =  pt.get <float>("search_space.driver_mutation_rate.true");
};



bool IsSettingFileAvailable(char * setting_file_name) {
  std::string s_setting_file_name(setting_file_name);
  //cout << check_etalon << endl;
  std::size_t pos = s_setting_file_name.find("model");
  //cout<<NUM<<endl;
  //cout << "Size of std::size_t : " << sizeof(std::size_t) << endl;
  //system("pause");
  //cout << "ind =" << ind << endl;
  if (pos != (std::string::npos)) {
        return true;
  } else {
        return false;
  };

};

void ThrowErrorArgumentNumber() {
  char err_mes[256] ="Mistake in the number of arguments:\n";
  char buffer[10];
  strcat(err_mes , buffer);
  strcat(err_mes ,"2 arguments for etalon probe or\n");
  strcat(err_mes ,"5 arguments for statistics\n");
  strcat(err_mes ,"6 arguments for the parameter inference simulation\n");
  strcat(err_mes ,"8 arguments for the distance comparison\n");
  err(err_mes);
};



void doesEtalonProbPart(BasicSimParameters & pars) {
	getEtalonPartParameters(pars);
	initSimGlobalValues(pars);
	if (run(pars, 1)) {
		//Initialize probe area (shoud be changed later)
		initMinMaxBordars(pars);
		//CreateInferenceConfFile(pars);
		SetEtalonProbe(pars);
	} else {
		std::cout << "num_of_attemps>"<<
		pars.threshold_simulation_num << std::endl;
	};
};

void doesInferencePart(AllSimParameters & pars) {
	getInferencePartParameters(pars.basic_sim_parameters,
														 pars.inference_parameters);
	initSimGlobalValues(pars.basic_sim_parameters); //TODO change this
	if (run(pars.basic_sim_parameters, 0)) {
		initMinMaxBordars(pars.basic_sim_parameters);
		DoesInferenceAnalyze(pars);
	} else {
		std::cout << "num_of_attemps>" <<
		pars.basic_sim_parameters.threshold_simulation_num << std::endl;
	};
};

void DoesDistInferencePart(AllSimParameters & pars) {
	initSimGlobalValues(pars.basic_sim_parameters); //TODO change this
	if (run(pars.basic_sim_parameters, 0)) {
		initMinMaxBordars(pars.basic_sim_parameters);
		DoesInferenceAnalyze(pars);
	} else {
		std::cout << "num_of_attemps>" <<
		pars.basic_sim_parameters.threshold_simulation_num << std::endl;
	};
};



void initParameterFileNamePaths(const char *argv[], BasicSimParameters & pars){
	pars.related_path_to_par_file = argv[1];
	pars.par_file_name = argv[2];
	std::string full_file_name = pars.related_path_to_par_file + pars.par_file_name;
	if (!fileExists(full_file_name.c_str())){
		std::cout << " file with the name "<<  full_file_name <<
                 " doesn't exist"<< std::endl;
		err("mistake in name settings");
	};
};

void initSearchSpacePar(const char *argv[], AllSimParameters & pars){
	 pars.basic_sim_parameters.evolution.driver_adv = atof(argv[3]);
	 pars.basic_sim_parameters.evolution.mutation_rate = atof(argv[4]);
	 pars.basic_sim_parameters.evolution.driver_mutation_rate = atof(argv[5]);
};

inline void etalonProbePart(const char *argv[]) {
	BasicSimParameters pars;
	initParameterFileNamePaths(argv, pars);
	doesEtalonProbPart(pars);
};

inline void simulationPart(const char *argv[]) {
	AllSimParameters pars;
	initParameterFileNamePaths(argv, pars.basic_sim_parameters);
	getAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
	initSearchSpacePar(argv, pars);
	pars.basic_sim_parameters.seed = atoi(argv[6]);
	doesInferencePart(pars);
};



void distComparisonPart(const char *argv[]) {
	AllSimParameters pars;
	initParameterFileNamePaths(argv, pars.basic_sim_parameters);
	getAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
// number of different search space points (general number 2 * point_number + 1)
	int point_number = atoi(argv[4]);
	int number_of_repeats = atoi(argv[5]); // for each search space point

	float left_scale = atof(argv[6]); // value lies in [0.0;1.0]
	float right_scale = atof(argv[7]); // value lies in [0.0;1.0]

	float d_a = pars.basic_sim_parameters.evolution.driver_adv;
	float m_r = pars.basic_sim_parameters.evolution.mutation_rate;
	float d_m_r = pars.basic_sim_parameters.evolution.driver_mutation_rate;
	switch(atoi(argv[8])) {
	case 0: {
	  			for(int i_dr_av = - point_number; i_dr_av <= point_number; i_dr_av++){
						if (i_dr_av < 0) pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + (d_a * left_scale * i_dr_av) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + ((1-d_a) * right_scale * i_dr_av) / (point_number * 1.0);
						for (int i = 0; i < number_of_repeats; i++){
							pars.basic_sim_parameters.seed++;
							DoesDistInferencePart(pars);
						};
	  			};
	  			pars.basic_sim_parameters.evolution.driver_adv = d_a;
	  			pars.basic_sim_parameters.evolution.mutation_rate = m_r;
	  			pars.basic_sim_parameters.evolution.driver_mutation_rate = d_m_r;
	  			for(int i_mu_r = - point_number; i_mu_r <= point_number; i_mu_r++){
						if (i_mu_r < 0) pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (m_r * left_scale * i_mu_r) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (( 1 - m_r) * right_scale * i_mu_r) / (point_number * 1.0);
						for (int i = 0; i < number_of_repeats; i++){
							 pars.basic_sim_parameters.seed++;
							 DoesDistInferencePart(pars);
						};
	  			};
					pars.basic_sim_parameters.evolution.driver_adv = d_a;
	  			pars.basic_sim_parameters.evolution.mutation_rate = m_r;
	  			pars.basic_sim_parameters.evolution.driver_mutation_rate = d_m_r;
	  			for(int i_dr_mu_r = - point_number; i_dr_mu_r <= point_number; i_dr_mu_r++){
						if (i_dr_mu_r < 0) pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (d_m_r * left_scale * i_dr_mu_r) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (( 1 - d_m_r) * right_scale * i_dr_mu_r) / (point_number * 1.0);
						for (int i = 0; i < number_of_repeats; i++){
							pars.basic_sim_parameters.seed++;
							DoesDistInferencePart(pars);
						};
					};
    } break;
	  case 1: {
		for(int i_dr_av = - point_number; i_dr_av <= point_number; i_dr_av++){
			if (i_dr_av < 0) pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + (d_a * left_scale * i_dr_av) / (point_number * 1.0);
			else pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + ((1-d_a) * right_scale * i_dr_av) / (point_number * 1.0);

			for(int i_mu_r = - point_number; i_mu_r <= point_number; i_mu_r++){
				if (i_mu_r < 0) pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (m_r * left_scale * i_mu_r) / (point_number * 1.0);
				else pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (( 1 - m_r) * right_scale * i_mu_r) / (point_number * 1.0);

				for(int i_dr_mu_r = - point_number; i_dr_mu_r <= point_number; i_dr_mu_r++){
					if (i_dr_mu_r < 0) pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (d_m_r * left_scale * i_dr_mu_r) / (point_number * 1.0);
					else pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (( 1 - d_m_r) * right_scale * i_dr_mu_r) / (point_number * 1.0);

					for (int i = 0; i < number_of_repeats; i++){
						pars.basic_sim_parameters.seed++;
						DoesDistInferencePart(pars);
					};
				};
			};
		};
	  } break;

	  case 2: {

	  			for(int i_dr_av = - point_number; i_dr_av <= point_number; i_dr_av++){
						if (i_dr_av < 0) pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + (d_a * left_scale * i_dr_av) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.driver_adv =
												 d_a + i_dr_av * 1.0;
						for (int i = 0; i < number_of_repeats; i++){
							pars.basic_sim_parameters.seed++;
							DoesDistInferencePart(pars);
						};
	  			};
	  			pars.basic_sim_parameters.evolution.driver_adv = d_a;
	  			pars.basic_sim_parameters.evolution.mutation_rate = m_r;
	  			pars.basic_sim_parameters.evolution.driver_mutation_rate = d_m_r;
	  			for(int i_mu_r = - point_number; i_mu_r <= point_number; i_mu_r++){
						if (i_mu_r < 0) pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (m_r * left_scale * i_mu_r) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.mutation_rate =
												 m_r + (( 1 - m_r) * right_scale * i_mu_r) / (point_number * 1.0);
						for (int i = 0; i < number_of_repeats; i++){
							pars.basic_sim_parameters.seed++;
							DoesDistInferencePart(pars);
						};
	  			};
					pars.basic_sim_parameters.evolution.driver_adv = d_a;
	  			pars.basic_sim_parameters.evolution.mutation_rate = m_r;
	  			pars.basic_sim_parameters.evolution.driver_mutation_rate = d_m_r;
	  			for(int i_dr_mu_r = - point_number; i_dr_mu_r <= point_number; i_dr_mu_r++){
						if (i_dr_mu_r < 0) pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (d_m_r * left_scale * i_dr_mu_r) / (point_number * 1.0);
						else pars.basic_sim_parameters.evolution.driver_mutation_rate =
												 d_m_r + (( 1 - d_m_r) * right_scale * i_dr_mu_r) / (point_number * 1.0);
						for (int i = 0; i < number_of_repeats; i++){
							pars.basic_sim_parameters.seed++;
							DoesDistInferencePart(pars);
						};
					};
    } break;
		default:{
			err("argument one or all");
    }
		};

};

inline void statisticsPart(const char *argv[]) {
	AllSimParameters pars;
	pars.basic_sim_parameters.related_path_to_par_file = argv[1];
	pars.basic_sim_parameters.par_file_name = argv[2];
	getAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
	std::string mode(argv[3]);
	if ( mode == "averaging") MakeAverageFile(pars);
};

inline void runProg(const int argc, const char *argv[]) {
// the number of arguments in command line and relevant regimes
  const int expected_arg_num_etalon_probe = 3;
  const int expected_arg_num_for_simulation = 7;
  const int expected_arg_num_for_dist_comparison = 9;
	const int expected_arg_num_for_statistics = 4;
// The number of parameters determine correspondent regime of the program
  switch(argc){
    case expected_arg_num_etalon_probe: {
			etalonProbePart(argv);
    } break;
    case expected_arg_num_for_simulation: {
 			simulationPart(argv);
    } break;
    case expected_arg_num_for_dist_comparison: {
		  distComparisonPart(argv);
		} break;
		case expected_arg_num_for_statistics: {
		} break;
    default:{
      ThrowErrorArgumentNumber();
    }
  }
}

int main(const int argc, const char *argv[]){
  runProg(argc, argv);
  return 0;
};


