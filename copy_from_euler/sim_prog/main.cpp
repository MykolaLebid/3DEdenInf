// Copyright 2017 Lebid Mykola all rights reserved

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <stdexcept>// out_of_range

using namespace std;
#include "params.h"
#include "classes.h"
#include "sample_processing.h"



#define __MAIN

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

//checking existence of a file
//#include <boost/filesystem.hpp>

//timer
#include <chrono>


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



void save_positions(char *name, float dz) {
  FILE *data=fopen(name,"w") ;
  for (unsigned int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    //Genotype *g=genotypes[cells[i].gen] ;
    if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<1e4)
        fprintf(data,"%d %d %d %d %d\n",i,
                int(cells[i].x+ll->r.x),
                int(cells[i].y+ll->r.y),
                int(cells[i].z+ll->r.z),
                cells[i].gen);
// fprintf(data,"%d %d %d %d %u\n",int(cells[i].x+ll->r.x),
// int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z),
// genotypes[cells[i].gen]->index) ;//WTF
  }
  fclose(data) ;
}

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

void SaveGenotypes(const char *name) {
  FILE *data=fopen(name,"w") ;
  for (unsigned int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d   %d   %d   (%d   %d)   %d\t"
              ,i,  g->prev_gen, g->num_snp, g->no_drivers,
               g->no_resistant, g->number) ;
      for (unsigned int j=0;j<g->sequence.size();j++) {
//mykola: we cut of additional bits for relevent number
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if ((g->sequence[j] & RESISTANT_PM)==0) {

                if ((g->sequence[j] & DRIVER_PM)==0) {
                    fprintf(data,"p %u  ",g->sequence[j]);
                } else {
                    fprintf(data,"d %u  ",(g->sequence[j] & (~DRIVER_PM)));
                };
            } else {
                fprintf(data,"r %u  ",(g->sequence[j] & (~RESISTANT_PM)));
            };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      };
      fprintf(data,"\n");
    };
  };

  fclose(data) ;
}



void save_most_abund_gens(char *name, unsigned int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (unsigned int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (unsigned int j=0;j<gg->sequence.size();j++) {
        if ((gg->sequence[j]&L_PM)==most_abund[0]) r=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[1]) g=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[2]) b=1 ;
      }
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index);
    }
  }
  fclose(data);
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

    if (s>0) printf("%d \n",s);
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
    pars.probe.population_borders.max_x = maxx;
    pars.probe.population_borders.max_y = maxy;
    pars.probe.population_borders.max_z = maxz;
    pars.probe.population_borders.min_x = minx;
    pars.probe.population_borders.min_y = miny;
    pars.probe.population_borders.min_z = minz;
};

void ThrowErrorArgumentNumber() {
  char err_mes[256] ="Mistake in the number of arguments:\n";
  char buffer[10];
  strcat(err_mes , buffer);
  strcat(err_mes ,"2 arguments for etalon probe or\n");
  strcat(err_mes ,"6 arguments for the parameter inference simulation\n");
  err(err_mes);
};

void InitSimGlobalValues(const BasicSimParameters & pars) {
    RAND = pars.seed;
    //Name a directory with etalon_mut_vec.dat (etalon mutation vector)
    std::string str_1 = pars.etalon_result_file_name;
    char *cstr = new char[str_1.length() + 1];
    strcpy(cstr, str_1.c_str());
    NUM = cstr;
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
        cout << "number of attempts to conduct the simulation is more then "<<
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
			InitSimGlobalValues(pars.basic_sim_parameters);
			if (run(pars.basic_sim_parameters, 0)) {
				InitPopulationBorders(pars.basic_sim_parameters);
				DoesInferenceAnalyze(pars);
			} else {
				cout << "number of attempts to conduct the simulation is more then " <<
				pars.basic_sim_parameters.threshold_simulation_num << endl;
			};
    } else {
			cout << "file with the name "<<  pars.basic_sim_parameters.par_file_name <<
					    " doesn't exist"<< endl;
    };
};

void RunProg(int argc, char *argv[]) {
  // the number of arguments in command line
  const int expected_arg_num_etalon_probe = 3;
  const int expected_arg_num_for_simulation = 7;
  switch(argc){
    case expected_arg_num_etalon_probe: { // 1 case
				BasicSimParameters pars;
				pars.related_path_to_par_file = argv[1];
				pars.par_file_name = argv[2];
				DoesEtalonProbPart(pars);
    }; break;
    case expected_arg_num_for_simulation: { // 2 case
        AllSimParameters pars;
        pars.basic_sim_parameters.related_path_to_par_file = argv[1];
        pars.basic_sim_parameters.par_file_name = argv[2];
        pars.basic_sim_parameters.evolution.driver_adv = atof(argv[3]);
        pars.basic_sim_parameters.evolution.mutation_rate = atof(argv[4]);
        pars.basic_sim_parameters.evolution.driver_mutation_rate = atof(argv[5]);
        pars.basic_sim_parameters.seed = atoi(argv[6]);
				DoesInferencePart(pars);
    }; break;
    default:
      ThrowErrorArgumentNumber(); break; // 3 case
  };
};

//void CreateInferenceConfFile(Parameters & pars) {
//  string full_way_to_config_file = pars.catalogue_file_name+"/"+
//                                   "model_inference_info.info";
//  boost::property_tree::ptree pt;
//
//  pt.add("search_space.driver_advantage.min" , pars.evolution.driver_adv/10.0 );
//  pt.add("search_space.driver_advantage.max" , pars.evolution.driver_adv*10.0);
//  pt.add("search_space.driver_advantage.true" , pars.evolution.driver_adv );
//
//  pt.add("search_space.mutation_rate.min" , pars.evolution.mutation_rate/10.0 );
//  pt.add("search_space.mutation_rate.max" , pars.evolution.mutation_rate*10  );
//  pt.add("search_space.mutation_rate.true" , pars.evolution.mutation_rate);
//
//  pt.add("search_space.driver_mutation_rate.min" ,
//         pars.evolution.driver_mutation_rate/10.0 );
//  pt.add("search_space.driver_mutation_rate.max" ,
//         pars.evolution.driver_mutation_rate*10);
//  pt.add("search_space.driver_mutation_rate.true" ,
//         pars.evolution.driver_mutation_rate);
//
//  pt.add("algorithm_parameters.final_epsilon" , 0.1 );
//  pt.add("algorithm_parameters.p_acc_min" , 0.00010 );
//  pt.add("algorithm_parameters.parameter_population_number" , 100 );
//  pt.add("algorithm_parameters.final_epsilon" , 100 );
//
//
//  pt.add("basic_settings.file_name" , pars.catalogue_file_name );
//  pt.add("basic_settings.seed" ,std::chrono::system_clock::
//         now().time_since_epoch().count());
//  pt.add("Parameters.ProbeComparisonParameters.comparison_file_name",
//         "file_tree_comparison.dat");
//  pt.add("Parameters.ProbeComparisonParameters.log_file_name",
//				  "logs.dat");
//  pt.add("Parameters.ProbeComparisonParameters.catalogue_with_e_mut_file",
//          pars.catalogue_file_name);
//  pt.add("Parameters.ProbeComparisonParameters.e_mut_file_name",
//          pars.file_name);
//  pt.add("Parameters.ProbeComparisonParameters.e_seed",
//          pars.seed);
//
//  pt.add("Parameters.ProbeComparisonParameters.Probe.cell_num",
//          pars.probe.cell_num);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.is_random",
//      pars.probe.is_random);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_x",
//          pars.probe.delta_area.del_x);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_y",
//          pars.probe.delta_area.del_y);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_z",
//          pars.probe.delta_area.del_z);
//
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_x",
//          pars.probe.population_borders.min_x);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_y",
//          pars.probe.population_borders.min_y);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_z",
//          pars.probe.population_borders.min_z);
//
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_x",
//          pars.probe.population_borders.max_x);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_y",
//          pars.probe.population_borders.max_y);
//  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_z",
//          pars.probe.population_borders.max_z);
//
//  pt.add("Parameters.ProbeComparisonParameters.Evolution.driver_adv",
//          pars.evolution.driver_adv);
//  pt.add("Parameters.ProbeComparisonParameters.Evolution.mutation_rate",
//          pars.evolution.mutation_rate);
//  pt.add("Parameters.ProbeComparisonParameters.Evolution.driver_mutation_rate",
//          pars.evolution.driver_mutation_rate);
//  pt.add("Parameters.ProbeComparisonParameters.Evolution.stop_time_cell_num",
//          pars.evolution.stop_time_cell_num);
//	pt.add("Parameters.ProbeComparisonParameters.dist_type", 1);
//  pt.add("Parameters.ProbeComparisonParameters.piece_num", 1);
//  pt.add("Parameters.ProbeComparisonParameters.max_dist_probe_trees", 1000);
//  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_cell_diff", 0.1);
//  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_mut_diff", 0.2);
//  pt.add("Parameters.threshold_simulation_num", pars.threshold_simulation_num);
//  boost::property_tree::info_parser::write_info(full_way_to_config_file, pt);
//};







int main(int argc, char *argv[]) {
  SetPreprocessorMethod();
  RunProg(argc, argv);
  return 0;
}
