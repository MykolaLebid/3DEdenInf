// Copyright 2017 Lebid Mykola all rights reserved
#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "params.h"
#include "classes.h"
#include "sample_processing.h"

#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL)
  #error no method defined!
#endif

#if (defined(VON_NEUMANN_NEIGHBOURHOOD) && defined(MOORE_NEIGHBOURHOOD))
  #error both VON_NEUMANN_NEIGHBOURHOOD and MOORE_NEIGHBOURHOOD defined!
#endif

#if (!defined(VON_NEUMANN_NEIGHBOURHOOD) && !defined(MOORE_NEIGHBOURHOOD))
  #error neither VON_NEUMANN_NEIGHBOURHOOD nor MOORE_NEIGHBOURHOOD defined!
#endif


bool run(const bool is_etalon_sim,
				 BasicSimParameters & pars,
				 SimData & sim_data)
{//init file system
	std::string path = pars.related_path_to_par_file + pars.etalon_result_catalogue;
	init(is_etalon_sim, path);
	//init data
	initSimData(pars.evolution.driver_adv,
							pars.evolution.driver_mutation_rate,
							sim_data);
	int s = 0;	// number of repeats
	// choose the regime (0 - well_mixed or 1 - 3D (spatial))
	switch (pars.evolution.simulation_regime) {
		case 0:{// well mixed
			while ((mainProgWellMixed(pars.evolution, sim_data)==1) &&
				(s < pars.threshold_simulation_num)) {
				initSimData(pars.evolution.driver_adv,
									  pars.evolution.driver_mutation_rate,
										sim_data);
				s++;
			} // initial growth until max size is reached
		} break;

		case 1:{// 3D
			while ((mainProg3D(pars.evolution, sim_data)==1) &&
				(s < pars.threshold_simulation_num)) {
					initSimData(pars.evolution.driver_adv, pars.evolution.driver_mutation_rate,
								sim_data);
					s++;
				} // initial growth until max size is reached
		} break;
		default:
			err("Parameter problem. Default case. pars.evolution.simulation_regime");
		break;
	}
	if (s < pars.threshold_simulation_num) return true;
	else return false;
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

void throwErrorArgumentNumber() {
  char err_mes[256] ="Mistake in the number of arguments:\n";
  char buffer[10];
  strcat(err_mes , buffer);
  strcat(err_mes ,"2 arguments for etalon probe or\n");
  strcat(err_mes ,"5 arguments for statistics\n");
  strcat(err_mes ,"6 arguments for the parameter inference simulation\n");
  strcat(err_mes ,"8 arguments for the distance comparison\n");
  err(err_mes);
};

void doesEtalonProbPart(BasicSimParameters & pars)
{
	getEtalonPartParameters(pars);
	SimData sim_data;
	bool is_etalon_prob = true;

	if (run(is_etalon_prob, pars, sim_data)) {
		setEtalonProbe(pars, sim_data);
	} else {
		std::cout << "num_of_attemps>"<< pars.threshold_simulation_num << std::endl;
	};
	cleanSimData(sim_data);
};

void doesInferencePart(AllSimParameters & pars)
{
	getInferencePartParameters(pars.basic_sim_parameters,
														 pars.inference_parameters);
	//initSimGlobalValues(pars.basic_sim_parameters); //TODO change this
	SimData sim_data;
	bool is_etalon_prob = false;
	if (run(is_etalon_prob, pars.basic_sim_parameters,  sim_data)) {
		//doesComparativeAnalysis(pars, sim_data);
	} else {
		std::cout << "num_of_attemps>" <<
		pars.basic_sim_parameters.threshold_simulation_num << std::endl;
	};
	cleanSimData(sim_data);
};

void doesDistInferencePart(AllSimParameters & pars)
{
	_srand48(pars.basic_sim_parameters.seed); //init seed
	SimData sim_data;
  bool is_etalon_prob = false;
	if (run(is_etalon_prob, pars.basic_sim_parameters, sim_data)) {
		//doesComparativeAnalysis(pars, sim_data);
	} else {
		std::cout << "num_of_attemps>" <<
		pars.basic_sim_parameters.threshold_simulation_num << std::endl;
	};
	cleanSimData(sim_data);
};

void initParameterFileNamePaths(const char *argv[], BasicSimParameters & pars)
{
	pars.related_path_to_par_file = argv[1];
	pars.par_file_name = argv[2];
	std::string full_file_name = pars.related_path_to_par_file +
															 pars.par_file_name;
	if (!fileExists(full_file_name.c_str())){
		std::cout << " file with the name "<< full_file_name <<
                 " doesn't exist"<< std::endl;
		err("mistake in name settings");
	};
};

void initSearchSpacePar(const char *argv[], AllSimParameters & pars)
{
	 pars.basic_sim_parameters.evolution.driver_adv = atof(argv[3]);
	 pars.basic_sim_parameters.evolution.mutation_rate = atof(argv[4]);
	 pars.basic_sim_parameters.evolution.driver_mutation_rate = atof(argv[5]);
};

inline void etalonProbePart(const char *argv[])
{
	BasicSimParameters pars;
	initParameterFileNamePaths(argv, pars);
	doesEtalonProbPart(pars);
};

inline void simulationPart(const char *argv[])
{
	AllSimParameters pars;
	initParameterFileNamePaths(argv, pars.basic_sim_parameters);
	getAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
	initSearchSpacePar(argv, pars);
	pars.basic_sim_parameters.seed = atoi(argv[6]);
	doesInferencePart(pars);
};

enum class ParameterEnum{
	driver_adv,
	mutation_rate,
	driver_mutation_rate
};

float getParValue(const ParameterEnum par_enum,
									AllSimParameters&   pars)
{
	float par_value;
	switch (par_enum) {
	case ParameterEnum::driver_adv:
		par_value = pars.basic_sim_parameters.evolution.driver_adv;
	break;
	case ParameterEnum::mutation_rate:
		par_value = pars.basic_sim_parameters.evolution.mutation_rate;
	break;
	case ParameterEnum::driver_mutation_rate:
		par_value = pars.basic_sim_parameters.evolution.driver_mutation_rate;
	break;
	default:
		err("Parameter problem. Default case");
	break;
	};
	return par_value;
};

void changeInitialParValue(const ParameterEnum par_enum,
												   const float 				 par_value,
												   AllSimParameters&   pars)
{
	switch (par_enum) {
	case ParameterEnum::driver_adv:
		pars.basic_sim_parameters.evolution.driver_adv           = par_value;
	break;
	case ParameterEnum::mutation_rate:
		pars.basic_sim_parameters.evolution.mutation_rate        = par_value;
	break;
	case ParameterEnum::driver_mutation_rate:
		pars.basic_sim_parameters.evolution.driver_mutation_rate = par_value;
	break;
	default:
		err("Parameter problem. Default case");
	break;
	};
};



void repeatDistInferencePart(const unsigned int number_of_repeats,
														 AllSimParameters&  pars)
{
	//#include <ctime>
	//clock_t start_time = clock(); //Start timer
	//double last;
	//double seconds_passed_from_start = 0;

	for (unsigned int i = 0; i < number_of_repeats; i++){
		pars.basic_sim_parameters.seed += 200;

		doesDistInferencePart(pars);
	};
	//last = seconds_passed_from_start;
	//	seconds_passed_from_start = (clock() - start_time);
	//	std::cout <<"seconds from last = " << (seconds_passed_from_start - last) <<'\n';
}

float returnConditionalParValue(const int 	index,
																const int 	point_number,
																const float par_value,
																const float left_scale,
																const float right_scale,
																const bool	scale_to_one
																)
{
	if (scale_to_one) {
		if (index < 0)
			return par_value + (par_value * left_scale * index) / (point_number * 1.0);
		else
			return par_value + ((1 - par_value) * right_scale * index) / (point_number * 1.0);
	} else {
		if (index < 0)
			return par_value + (par_value * left_scale * index) / (point_number * 1.0);
		else
			return par_value + (par_value * right_scale * index) / (point_number * 1.0);
	};


}

void runSimulationLoopOnePar(const ParameterEnum	par_enum,
													   const unsigned int		point_number,
													   const unsigned int		repeat_number,
													   const float					left_scale,
													   const float					right_scale,
													   AllSimParameters&    pars)
{
	bool scale_to_one = pars.inference_parameters.comparison_parameters.par_scale_to_one;
	float par_value 	= getParValue(par_enum, pars);

	for(int i = - point_number ; i <= (int) point_number; i++) {
		float new_par_value = returnConditionalParValue(i, point_number, par_value,
																										left_scale, right_scale, scale_to_one);
		changeInitialParValue(par_enum, new_par_value, pars);
		repeatDistInferencePart(repeat_number, pars);
		std::cout<<"simulation " << i + point_number + 1 << " is done" <<"\n";
	};
}


void runSimulationLoopAllPar(const unsigned int	 	point_num,
													   const unsigned int	 	repeat_number,
													   const float          left_scale,
													   const float					right_scale,
													   AllSimParameters&    pars)
{
	bool scale_to_one = pars.inference_parameters.comparison_parameters.par_scale_to_one;
	float d_a   = pars.basic_sim_parameters.evolution.driver_adv;
	float m_r   = pars.basic_sim_parameters.evolution.mutation_rate;
  float d_m_r = pars.basic_sim_parameters.evolution.driver_mutation_rate;

	for(int i_dr_av = - point_num; i_dr_av <= (int) point_num; i_dr_av++){
		pars.basic_sim_parameters.evolution.driver_adv =
				 returnConditionalParValue(i_dr_av, point_num, d_a,
																	 left_scale, right_scale, scale_to_one);
		for(int i_mu_r = - point_num; i_mu_r <= (int) point_num; i_mu_r++){
			pars.basic_sim_parameters.evolution.mutation_rate =
					returnConditionalParValue(i_mu_r, point_num, m_r,
																		left_scale, right_scale, scale_to_one);
			for(int i_dr_mu_r = - point_num; i_dr_mu_r <= (int) point_num; i_dr_mu_r++){
				pars.basic_sim_parameters.evolution.driver_mutation_rate =
						returnConditionalParValue(i_dr_mu_r, point_num, d_m_r,
																			left_scale, right_scale, scale_to_one);
						repeatDistInferencePart(repeat_number, pars);
			}
		}
	}

}

void distComparisonPart(const char *argv[])
{
	AllSimParameters pars;
	initParameterFileNamePaths(argv, pars.basic_sim_parameters);
	getAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
// number of different search space points (general number 2 * point_number + 1)
  unsigned int point_number 			= atoi(argv[4]); // general number
	unsigned int number_of_repeats	= atoi(argv[5]); // for each search space point
	float left_scale								= atof(argv[6]); // value lies in [0.0;1.0]
	float right_scale  							= atof(argv[7]); // value lies in [0.0;1.0]
	unsigned int mode 							= atoi(argv[8]); // 1 or 2
	switch(mode) {
		case 0: {
			pars.inference_parameters.comparison_parameters.par_scale_to_one = true;
//			runSimulationLoopOnePar(ParameterEnum::driver_adv, point_number,
//															number_of_repeats, left_scale, right_scale, pars);
	  	runSimulationLoopOnePar(ParameterEnum::mutation_rate, point_number,
															number_of_repeats, left_scale, right_scale, pars);
//	  	runSimulationLoopOnePar(ParameterEnum::driver_mutation_rate, point_number,
//															number_of_repeats, left_scale, right_scale, pars);
		} break;
	  case 1: {
			pars.inference_parameters.comparison_parameters.par_scale_to_one = true;
			runSimulationLoopAllPar(point_number, number_of_repeats,
															left_scale, right_scale, pars);
	  } break;
	  case 2: {
			pars.inference_parameters.comparison_parameters.par_scale_to_one = false;
//	  	runSimulationLoopOnePar(ParameterEnum::driver_adv, point_number,
//															number_of_repeats, left_scale, right_scale, pars);
	  	runSimulationLoopOnePar(ParameterEnum::mutation_rate, point_number,
															number_of_repeats, left_scale, right_scale, pars);
//	  	runSimulationLoopOnePar(ParameterEnum::driver_mutation_rate, point_number,
//															number_of_repeats, left_scale, right_scale, pars);
    } break;
		default:{
			err("argument one or all");
    };
	};
}

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
			statisticsPart(argv);
		} break;
    default:{
      throwErrorArgumentNumber();
    }
  }
}

int main(const int argc, const char *argv[]){
  runProg(argc, argv);
  return 0;
};
