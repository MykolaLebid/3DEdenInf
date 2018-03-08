// Copyright 2017 Lebid Mykola all rights reserved

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
//#include <fstream>
#include <stdexcept>// out_of_range

//using namespace std;
#include "params.h"
#include "classes.h"
#include "sample_processing.h"


#define __MAIN

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

//checking existence of a file
//#include <boost/filesystem.hpp>

//timer
//#include <chrono>


#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL)
  #error no method defined!
#endif

#if (defined(VON_NEUMANN_NEIGHBOURHOOD) && defined(MOORE_NEIGHBOURHOOD))
  #error both VON_NEUMANN_NEIGHBOURHOOD and MOORE_NEIGHBOURHOOD defined!
#endif

#if (!defined(VON_NEUMANN_NEIGHBOURHOOD) && !defined(MOORE_NEIGHBOURHOOD))
  #error neither VON_NEUMANN_NEIGHBOURHOOD nor MOORE_NEIGHBOURHOOD defined!
#endif

extern char *NUM;
extern const char *NUM_E;
extern std::string STR;
extern int RAND, sample, treatment, max_size, nsam;
extern double tt ;//time

int max_size;
float driver_adv;
float driver_prob;
float gama;
float death0, growth0=1;

int num_adv;

int sample=0;



int minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;

//search_min_max
void SearchMinMax(){
//Mykola search of min and max for all cells in all lesions
  for (unsigned int i=0; i < cells.size(); i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  };
  maxx++ ; maxy++;
}






bool run(BasicSimParameters & pars, bool is_etalon_sim) {
   std::string path = pars.related_path_to_par_file +
											pars.etalon_result_catalogue;
	// std::cout << path << std::endl;
   init(is_etalon_sim, path);
   int s = 0;
   reset();
   while ((main_proc(pars.evolution.stop_time_cell_num,
                     2,-1, -1)==1) &&
          (s < pars.threshold_simulation_num)) {
        s++ ; reset();
    } // initial growth until max size is reached

    //if (s>0) printf("%d \n",s);
    SearchMinMax();
    if (s < pars.threshold_simulation_num) return true;
    else return false;

};

void init_min_max()
{
    minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;
}

Borders GetPopulationBorders() {
  return  {minx, maxx,
           miny, maxy,
           minz, maxz};
};


void InitParamSearch (char *name) {
    //read parameters from file "model_info.info"
    std::string config = string(name);
    //cout<<config<<endl;

    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(config, pt);

    RAND = pt.get<int> ("basic_settings.seed");
    //Name a directory with etalon_mut_vec.dat (etalon mutation vector)
    std::string str_1 = pt.get<string> ("basic_settings.file_name");
    char *cstr = new char[str_1.length() + 1];
    strcpy(cstr, str_1.c_str());
    NUM = cstr;
    // do stuff
    //cout<< "NUM" << NUM << endl;
    //cout<< "str.c_str()="<<str.c_str()<<endl;
    //cout<< "(char*) str.c_str()="<<(char*) str.c_str()<<endl;
    max_size    =  pt.get <float>("evolution_settings.stopping_time");
    driver_adv  =  pt.get <float>("search_space.driver_advantage.true");
    gama        =  pt.get <float>("search_space.mutation_rate.true");
    driver_prob =  pt.get <float>("search_space.driver_mutation_rate.true");
};


void SetPreprocessorMethod() {
#if defined(GILLESPIE)
 // cout <<"method: GILLESPIE\n" ;
#endif
#if defined(FASTER_KMC)
 // cout <<"method: FASTER_KMC\n" ;
#endif
#if defined(NORMAL)
 // cout <<"method: NORMAL\n" ;
#endif
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

void InitPopulationBorders(BasicSimParameters & pars) {
    pars.probes.population_borders.max_x = maxx;
    pars.probes.population_borders.max_y = maxy;
    pars.probes.population_borders.max_z = maxz;
    pars.probes.population_borders.min_x = minx;
    pars.probes.population_borders.min_y = miny;
    pars.probes.population_borders.min_z = minz;
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

void InitSimGlobalValues(const BasicSimParameters & pars) {
    RAND = pars.seed;
    //Name a directory with etalon_mut_vec.dat (etalon mutation vector)
    //std::string str_1 = pars.etalon_result_file_name;
    //std::string str_1 = pars.etalon_result_file_name;
    //char *cstr = new char[str_1.length() + 1];
    //strcpy(cstr, str_1.c_str());
    //NUM = cstr;
    max_size = pars.evolution.stop_time_cell_num;
    driver_adv  =  pars.evolution.driver_adv;
    gama  =  pars.evolution.mutation_rate;
    driver_prob =  pars.evolution.driver_mutation_rate;
		growth0 = 1;
    death0 = growth0 * (1 - pars.evolution.driver_adv);
    _srand48(pars.seed);
};


void DoesEtalonProbPart(BasicSimParameters & pars) {
		std::string full_file_name = pars.related_path_to_par_file +
																 pars.par_file_name;
		if (FileExists(full_file_name.c_str())){
      GetEtalonPartParameters(pars);
      InitSimGlobalValues(pars);
      if (run(pars, 1)) {
        InitPopulationBorders(pars);
        //Initialize probe area (shoud be changed later)
        SetEtalonProbe(pars);
        //CreateInferenceConfFile(pars);
      } else {
        cout << "num_of_attemps>"<<
        pars.threshold_simulation_num << endl;
      };
    } else {
      cout << "file with the name "<<  pars.par_file_name <<
              " doesn't exist"<< endl;
       //static_cast<const void*>()
    };
};

void DoesInferencePart(AllSimParameters & pars) {
		std::string full_file_name = pars.basic_sim_parameters.related_path_to_par_file
																 + pars.basic_sim_parameters.par_file_name;
    if (FileExists(full_file_name.c_str())){
      GetInferencePartParameters(pars.basic_sim_parameters,
																 pars.inference_parameters);
			InitSimGlobalValues(pars.basic_sim_parameters); //TODO change this
			if (run(pars.basic_sim_parameters, 0)) {
				InitPopulationBorders(pars.basic_sim_parameters);
				DoesInferenceAnalyze(pars);
			} else {
				cout << "num_of_attemps>" <<
				pars.basic_sim_parameters.threshold_simulation_num << endl;
			};
    } else {
			cout << "file with the name "<<  pars.basic_sim_parameters.par_file_name <<
					    " doesn't exist"<< endl;
    };
};

void DoesDistInferencePart(AllSimParameters & pars) {
			InitSimGlobalValues(pars.basic_sim_parameters); //TODO change this
			if (run(pars.basic_sim_parameters, 0)) {
				InitPopulationBorders(pars.basic_sim_parameters);
				DoesInferenceAnalyze(pars);
			} else {
				cout << "num_of_attemps>" <<
				pars.basic_sim_parameters.threshold_simulation_num << endl;
			};
};

void InitSearchSpacePar(char *argv[], AllSimParameters & pars){
	 pars.basic_sim_parameters.evolution.driver_adv = atof(argv[3]);
	 pars.basic_sim_parameters.evolution.mutation_rate = atof(argv[4]);
	 pars.basic_sim_parameters.evolution.driver_mutation_rate = atof(argv[5]);
};

void DistComparison(int argc, char *argv[]){
	AllSimParameters pars;
	pars.basic_sim_parameters.related_path_to_par_file = argv[1];
	pars.basic_sim_parameters.par_file_name = argv[2];
	std::string full_file_name = pars.basic_sim_parameters.related_path_to_par_file
															 + pars.basic_sim_parameters.par_file_name;
	if (FileExists(full_file_name.c_str())){
		GetAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
		int point_number = atoi(argv[4]); // number of different search space points (general number 2 * point_number + 1)
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
};

void RunProg(int argc, char *argv[]) {
  // the number of arguments in command line
  const int expected_arg_num_etalon_probe = 3;
  const int expected_arg_num_for_statistics = 4;
  const int expected_arg_num_for_simulation = 7;
  const int expected_arg_num_for_dist_comparison = 9;
  switch(argc){
    case expected_arg_num_etalon_probe: {
				BasicSimParameters pars;
				pars.related_path_to_par_file = argv[1];
				pars.par_file_name = argv[2];
				DoesEtalonProbPart(pars);
    } break;
    case expected_arg_num_for_simulation: {
        AllSimParameters pars;
				pars.basic_sim_parameters.related_path_to_par_file = argv[1];
	      pars.basic_sim_parameters.par_file_name = argv[2];
	      GetAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
	      InitSearchSpacePar(argv, pars);
	      pars.basic_sim_parameters.seed = atoi(argv[6]);
				DoesInferencePart(pars);
    } break;
    case expected_arg_num_for_dist_comparison: {
		  DistComparison(argc, argv);
		} break;
		case expected_arg_num_for_statistics: {
        AllSimParameters pars;
				pars.basic_sim_parameters.related_path_to_par_file = argv[1];
	      pars.basic_sim_parameters.par_file_name = argv[2];
	      GetAllParameters(pars.basic_sim_parameters, pars.inference_parameters);
        std::string mode(argv[3]);
				if ( mode == "averaging") MakeAverageFile(pars);

		} break;
    default:{
      ThrowErrorArgumentNumber();
    }
  }
}



int main(int argc, char *argv[]){
  SetPreprocessorMethod();
  RunProg(argc, argv);
  return 0;
};
