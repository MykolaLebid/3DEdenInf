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
  fclose(data) ;

}

bool run(Parameters pars) {
   init(pars.is_simulation_for_etalon_probe,
        pars.catalog_file_name);
   int s = 0 ;
   reset();
   while ((main_proc(pars.evo_pars.stop_time_cell_num,
                     2,-1, -1)==1) &&
          (s<pars.threshold_simulation_num)) {
        s++ ; reset();
    } // initial growth until max size is reached

    if (s>0) printf("resetted %d times\n",s);
    SearchMinMax();
    if (s<pars.threshold_simulation_num) return true;
    else return false;
//#endif
}

void init_min_max()
{
    minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;
}

Borders GetPopulationBorders() {
  return  {minx, maxx,
           miny, maxy,
           minz, maxz};
};


//initializes parameters for for etalon probe
//
void GetEtalonProbePar(Parameters & pars) {
    //read parameters from file "model_info.info"
    //std::string config = string(name);
    //cout<<config<<endl;
    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(pars.par_file_name, pt);
// delete later
///////////////////////////////////////////////////////////////////////////////
    RAND = pt.get<int> ("Parameters.seed");                                  //
    //Name a directory with etalon_mut_vec.dat (etalon mutation vector)      //
    std::string str_1 = pt.get<string> ("Parameters.file_name");             //
    char *cstr = new char[str_1.length() + 1];                               //
    strcpy(cstr, str_1.c_str());                                             //
    NUM = cstr;
    // do stuff
    //cout<< "NUM" << NUM << endl;
    //cout<< "str.c_str()="<<str.c_str()<<endl;
    //cout<< "(char*) str.c_str()="<<(char*) str.c_str()<<endl;
    max_size    =  pt.get <int> ("Parameters.Evolution.stop_time_cell_num");
    driver_adv  =  pt.get <float>                                             //
    ("Parameters.Evolution.driver_adv");                                      //
    gama        =  pt.get <float>                                             //
    ("Parameters.Evolution.mutation_rate");                                   //
    driver_prob =  pt.get <float>                                             //
    ("Parameters.Evolution.driver_mutation_rate");                            //
////////////////////////////////////////////////////////////////////////////////
    pars.seed = pt.get<int> ("Parameters.seed");
    pars.file_name = pt.get<string> ("Parameters.file_name");
    pars.catalog_file_name = pt.get<string> ("Parameters.catalog_file_name");
    pars.is_simulation_for_etalon_probe = true;
    pars.threshold_simulation_num =
      pt.get<int> ("Parameters.threshold_simulation_num");
    pars.evo_pars.driver_adv =
      pt.get <float> ("Parameters.Evolution.driver_adv");
    pars.evo_pars.mutation_rate =
      pt.get <float>("Parameters.Evolution.mutation_rate");
    pars.evo_pars.driver_mutation_rate =
      pt.get <float>("Parameters.Evolution.driver_mutation_rate");
    pars.evo_pars.stop_time_cell_num =
      pt.get <int> ("Parameters.Evolution.stop_time_cell_num");
    pars.probe_pars.cell_num =
      pt.get <int> ("Parameters.Probe.cell_num");
    pars.probe_pars.is_random =
			pt.get <bool> ("Parameters.Probe.is_random");
    pars.probe_pars.del_area.del_x =
      pt.get <int> ("Parameters.Probe.DeltaArea.del_x");
    pars.probe_pars.del_area.del_y =
      pt.get <int> ("Parameters.Probe.DeltaArea.del_y");
    pars.probe_pars.del_area.del_z =
      pt.get <int> ("Parameters.Probe.DeltaArea.del_z");
    growth0 = 1;

    death0 = growth0 * (1 - pars.evo_pars.driver_adv);
    _srand48(pars.seed);

};

//initializes parameters for for inference part
//
void GetInferencePar(Parameters & pars) {
  //read parameters from file "model_info.info"
  //std::string config = string(name);
  //cout<<config<<endl;

  boost::property_tree::ptree pt;
  boost::property_tree::info_parser::read_info(pars.par_file_name, pt);
  // begin* Name a directory with etalon_mut_vec.dat (etalon mutation vector)
  pars.catalog_file_name = pt.get<string>
		("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file");
  char *cstr = new char[pars.catalog_file_name.length() + 1];
  strcpy(cstr, pars.catalog_file_name.c_str());
  // end*
  pars.evo_pars.stop_time_cell_num = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Evolution.stop_time_cell_num");


  // delete later
  ///////////////////////////////////////////////////////////////////////////
  driver_adv  =  pars.evo_pars.driver_adv;                          //
  gama        =  pars.evo_pars.mutation_rate;                       //
  driver_prob =  pars.evo_pars.driver_mutation_rate;                //
  RAND        =  pars.seed ;
  NUM         =  cstr;                                                     //
  max_size    =  pars.evo_pars.stop_time_cell_num;
  ///////////////////////////////////////////////////////////////////////////

  pars.file_name = pt.get<string>
    ("Parameters.ProbeComparisonParameters.e_mut_file_name");
  pars.etalon.comparison_file_name = pt.get<string>
    ("Parameters.ProbeComparisonParameters.comparison_file_name");
  pars.etalon.cataloge_with_e_mut_file = pt.get<string>
    ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file");
  pars.etalon.e_mut_file_name = pt.get<string>
    ("Parameters.ProbeComparisonParameters.e_mut_file_name");

  pars.probe_pars.cell_num = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.cell_num");
	pars.probe_pars.is_random = pt.get <bool>
		("Parameters.ProbeComparisonParameters.Probe.is_random");
  pars.probe_pars.del_area.del_x = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_x");
  pars.probe_pars.del_area.del_y = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_y");
  pars.probe_pars.del_area.del_z = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_z");
  pars.etalon.piece_num = pt.get<int>
    ("Parameters.ProbeComparisonParameters.piece_num");
  pars.etalon.dist_type = pt.get<unsigned int>
    ("Parameters.ProbeComparisonParameters.dist_type");
  pars.etalon.max_dist_probe_trees = pt.get<float>
    ("Parameters.ProbeComparisonParameters.max_dist_probe_trees");
  pars.etalon.threshold_frac_cell_diff = pt.get<float>
    ("Parameters.ProbeComparisonParameters.threshold_frac_cell_diff");
  pars.etalon.threshold_frac_mut_diff = pt.get<float>
    ("Parameters.ProbeComparisonParameters.threshold_frac_mut_diff");
  pars.threshold_simulation_num = pt.get <int>
    ("Parameters.threshold_simulation_num");

    growth0 = 1;
    death0 = growth0 * (1 - pars.evo_pars.driver_adv);
    _srand48(pars.seed);
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

bool FileExists(const char *name) {
    ifstream my_file(name);
    if(my_file.fail()){ //File does not exist
        return false;
    } else { //otherwise, file exists
        return true;
    };
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

void InitPopulationBorders(Parameters & pars) {
    pars.probe_pars.population_borders.max_x = maxx;
    pars.probe_pars.population_borders.max_y = maxy;
    pars.probe_pars.population_borders.max_z = maxz;
    pars.probe_pars.population_borders.min_x = minx;
    pars.probe_pars.population_borders.min_y = miny;
    pars.probe_pars.population_borders.min_z = minz;
};

bool IsSimulationForEtalonProbe(char * setting_file_name) {
  std::string check_etalon(setting_file_name);
//
//  #if defined __linux
//    std::ifstream infile("./model_info.info");
//  #else
//     std::ifstream infile("model_info.info");
//  #endif
//  if (infile.good()){
  if (check_etalon.find("model_etalon_info") != std::string::npos) {
        return true;
  } else {
        return false;
  };
};
bool IsSimulationForInference(char * setting_file_name) {
  std::string check_etalon(setting_file_name);
  if (check_etalon.find("model_inference_info") != std::string::npos) {
        return true;
  } else {
        return false;
  };
};


void ThrowErrorArgumentNumber() {
  char err_mes[256] ="Mistake in the number of arguments:\n";
  char buffer[10];
  strcat(err_mes , buffer);
  strcat(err_mes ,"1 argument for etalon probe or\n");
  strcat(err_mes ,"5 arguments for the parameter inference simulation\n");
  err(err_mes);
};

Parameters InitParameters(int argc, char *argv[]) {
  Parameters pars;
  // the number of arguments in command line
  const int expected_arg_num_etalon_probe = 2;
  const int expected_arg_num_for_simulation = 6;
  switch( argc ){
    case expected_arg_num_etalon_probe: { // 1
      if (IsSimulationForEtalonProbe(argv[1])) {
          pars.is_simulation_for_etalon_probe = true;
          pars.par_file_name = argv[1];
          return pars;
      } else {
        err("MISTAKE Name format (etalon): model_etalon_info*.info\n");
      };
    }; break;
    case expected_arg_num_for_simulation: { // 2
      if (IsSimulationForInference(argv[1])) {
        pars.par_file_name = argv[1];
        pars.is_simulation_for_etalon_probe = false;
        pars.evo_pars.driver_adv = atof(argv[2]);
        pars.evo_pars.mutation_rate = atof(argv[3]);
        pars.evo_pars.driver_mutation_rate = atof(argv[4]);
        pars.seed = atof(argv[5]);
        return pars;
      } else {
        err("MISTAKE Name format (inference): model_inference_info*.info\n");
      };
    }; break;
    default:
      ThrowErrorArgumentNumber(); break; // 3
  };
  return pars;
};

void CreateInferenceConfFile(Parameters & pars) {
  string full_way_to_config_file = pars.catalog_file_name+"/"+
                                   "model_inference_info.info";
  boost::property_tree::ptree pt;

  pt.add("search_space.driver_advantage.min" , pars.evo_pars.driver_adv/10.0 );
  pt.add("search_space.driver_advantage.max" , pars.evo_pars.driver_adv*10.0);
  pt.add("search_space.driver_advantage.true" , pars.evo_pars.driver_adv );

  pt.add("search_space.mutation_rate.min" , pars.evo_pars.mutation_rate/10.0 );
  pt.add("search_space.mutation_rate.max" , pars.evo_pars.mutation_rate*10  );
  pt.add("search_space.mutation_rate.true" , pars.evo_pars.mutation_rate);

  pt.add("search_space.driver_mutation_rate.min" ,
         pars.evo_pars.driver_mutation_rate/10.0 );
  pt.add("search_space.driver_mutation_rate.max" ,
         pars.evo_pars.driver_mutation_rate*10);
  pt.add("search_space.driver_mutation_rate.true" ,
         pars.evo_pars.driver_mutation_rate);

  pt.add("algorithm_parameters.final_epsilon" , 0.1 );
  pt.add("algorithm_parameters.p_acc_min" , 0.00010 );
  pt.add("algorithm_parameters.parameter_population_number" , 100 );


  pt.add("basic_settings.file_name" , pars.catalog_file_name );
  pt.add("basic_settings.seed" ,std::chrono::system_clock::
         now().time_since_epoch().count());
  pt.add("Parameters.ProbeComparisonParameters.comparison_file_name",
         "file_tree_comparison.dat");

  pt.add("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file",
          pars.catalog_file_name);
  pt.add("Parameters.ProbeComparisonParameters.e_mut_file_name",
          pars.file_name);
  pt.add("Parameters.ProbeComparisonParameters.e_seed",
          pars.seed);

  pt.add("Parameters.ProbeComparisonParameters.Probe.cell_num",
          pars.probe_pars.cell_num);
  pt.add("Parameters.ProbeComparisonParameters.Probe.is_random",
      pars.probe_pars.is_random);
  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_x",
          pars.probe_pars.del_area.del_x);
  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_y",
          pars.probe_pars.del_area.del_y);
  pt.add("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_z",
          pars.probe_pars.del_area.del_z);

  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_x",
          pars.probe_pars.population_borders.min_x);
  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_y",
          pars.probe_pars.population_borders.min_y);
  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.min_z",
          pars.probe_pars.population_borders.min_z);

  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_x",
          pars.probe_pars.population_borders.max_x);
  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_y",
          pars.probe_pars.population_borders.max_y);
  pt.add("Parameters.ProbeComparisonParameters.Probe.Area.max_z",
          pars.probe_pars.population_borders.max_z);

  pt.add("Parameters.ProbeComparisonParameters.Evolution.driver_adv",
          pars.evo_pars.driver_adv);
  pt.add("Parameters.ProbeComparisonParameters.Evolution.mutation_rate",
          pars.evo_pars.mutation_rate);
  pt.add("Parameters.ProbeComparisonParameters.Evolution.driver_mutation_rate",
          pars.evo_pars.driver_mutation_rate);
  pt.add("Parameters.ProbeComparisonParameters.Evolution.stop_time_cell_num",
          pars.evo_pars.stop_time_cell_num);
	pt.add("Parameters.ProbeComparisonParameters.dist_type", 1);
  pt.add("Parameters.ProbeComparisonParameters.piece_num", 1);
  pt.add("Parameters.ProbeComparisonParameters.max_dist_probe_trees", 1000);
  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_cell_diff", 0.1);
  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_mut_diff", 0.2);
  pt.add("Parameters.threshold_simulation_num", pars.threshold_simulation_num);
  boost::property_tree::info_parser::write_info(full_way_to_config_file, pt);
};

void DoesEtalonProbPart(Parameters & pars) {
    if (FileExists(pars.par_file_name.c_str())){
      GetEtalonProbePar(pars); // EtalonProbe
      if (run(pars)) {
        InitPopulationBorders(pars);
        //Initialize probe area (shoud be changed later)
        SetEtalonProbe(pars);
        CreateInferenceConfFile(pars);
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
void DoesInferencePart(Parameters & pars) {
    if (FileExists(pars.par_file_name.c_str())){
      GetInferencePar(pars);
        if (run(pars)) {

          InitPopulationBorders(pars);
          DoesInferenceAnalyze(pars);

        } else {
          cout << "number of attempts to conduct the simulation is more then "<<
          pars.threshold_simulation_num << endl;
        };
    } else {
      cout << "file with the name "<<  pars.par_file_name <<
              " doesn't exist"<< endl;
    };
};

int main(int argc, char *argv[]) {
  SetPreprocessorMethod();
  Parameters pars = InitParameters(argc, argv);
  if (pars.is_simulation_for_etalon_probe) {
    DoesEtalonProbPart(pars);
  } else {
    DoesInferencePart(pars);
  };
  return 0;
}
