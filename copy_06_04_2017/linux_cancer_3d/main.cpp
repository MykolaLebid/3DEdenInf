// Copyright 2017 Lebid Mykola all rights reserved

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <stdexcept>// out_of_range

using namespace std;
#include "classes.h"
#define __MAIN
#include "params.h"
#include "sample_processing.h"



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


float driver_adv;
float driver_prob;
float gama;
float death0, growth0=1;

int num_adv;

int sample=0;



int minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;



void save_positions(char *name, float dz) {
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
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

float save_2d_image(char *name, vecd li)// li direction of light
{
  int i,j,k;

  //Mykola search

  float density=1.*cells.size()/(float(maxx-minx)*float(maxy-miny)*float(maxz-minz)) ;
  if (cells.size()<1e3) return density ;


  float diam=pow(float(maxx-minx)*float(maxy-miny)*float(maxz-minz),1./3) ;
  //printf("density=%f\n",density) ;

  int nnn=(maxx-minx)*(maxy-miny) ;
  if (float(maxx-minx)*float(maxy-miny)>2e9) err("(maxx-minx)*(maxy-miny) too large",float(maxx-minx)*float(maxy-miny)) ;
  int *types=new int[nnn] ; // 2d array of cells' types (foremost ones only)
  short int *zbuf=new short int[nnn] ; // 2d array of cells' z positions (foremost ones only)
  BYTE *br=new BYTE[nnn] ; // 2d array of cells' brightness (foremost ones only)
  Sites **bit=new Sites*[nnn] ; // 3d array of bits: cell empty/occupied
  for (i=0;i<nnn;i++) { zbuf[i]=minz ; types[i]=-1 ; br[i]=255 ; bit[i]=new Sites(maxz-minz) ; }
  //printf("%d %d\t %d %d\n",minx,maxx,miny,maxy) ;
  //printf("%d x %d = %d\n",maxx-minx,maxy-miny,nnn) ;

  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int z=int(cells[i].z+ll->r.z) ;
    int adr=int(cells[i].y+ll->r.y-miny)*(maxx-minx)+int(cells[i].x+ll->r.x-minx) ;
    if (adr<0 || adr>=nnn) err("adr",adr) ;
    bit[adr]->set(z-minz) ;
    if (z>zbuf[adr]) { zbuf[adr]=z ; types[adr]=genotypes[cells[i].gen]->index ; }
  }

  normalize(li) ;
  float dmul=0.93 ;
  float range=-(0.916291/(density*log(dmul))) ;
  if (range>diam) { range=diam ; dmul=pow(0.4,1/(density*range)) ; }
  for (i=0;i<maxy-miny;i++)
    for (j=0;j<maxx-minx;j++) {
      float d=1 ;
      int adr=i*(maxx-minx)+j ;
      k=zbuf[adr]-minz ;
      for (float o=1;o<range;o++) {
        int il=int(i-li.y*o) ; int jl=int(j-li.x*o) ; int kl=int(k-li.z*o) ;
        if (il>0 && jl>0 && il<maxy-miny && jl<maxx-minx && kl>0 && kl<maxz-minz) {
          if (bit[il*(maxx-minx)+jl]->is_set(kl)) d*=dmul ;
        } else break ;
        if (d<0.4) { d=0.4 ; break ; }
      }
      br[adr]*=d ;
    }

  FILE *data=fopen(name,"w") ;
  for (i=0;i<maxy-miny;i++) {
    for (j=0;j<maxx-minx;j++) fprintf(data,"%d %d ",types[i*(maxx-minx)+j],br[i*(maxx-minx)+j]) ;
    fprintf(data,"\n") ;
  }
  fclose(data) ;
  delete [] types ; delete [] zbuf ; delete [] br ; delete [] bit ;
  return density ;
}


void SaveGenotypes(char *name) {
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d   %d   %d   (%d   %d)   %d\t"
              ,i,  g->prev_gen, g->num_snp, g->no_drivers,
               g->no_resistant, g->number) ;
      for (int j=0;j<g->sequence.size();j++) {
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



void save_most_abund_gens(char *name, int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (int j=0;j<gg->sequence.size();j++) {
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

/*
    fflush(stdout) ;

    save_data() ; //Mykola save_data in test_1

    sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
    FILE *er=fopen(name,"a") ;
    fprintf(er,"%d\t%d %lf\n",sample,s,tt) ;//Mykola save num of samples and time
    fclose(er) ;
    // save some more data

    int *snp_no=new int[L], *snp_drivers=new int[L] ; // array of SNPs abundances
    for (int i=0;i<L;i++) { snp_no[i]=snp_drivers[i]=0 ; }
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0)
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
          if (((genotypes[i]->sequence[j])&DRIVER_PM)) snp_drivers[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
        }
    }

    save_spatial(snp_no) ;
    //printf("time %lf\n",tt);
    //printf("saving PMs...\n") ;
    int most_abund[100] ;
    sprintf(name,"%s/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,max_size,0,most_abund) ;
*/
/*
    if (driver_adv>0 || driver_migr_adv>0) { //printf("saving driver PMs...\n") ; sprintf(name,"%s/drv_PMs_%d_%d.dat",NUM,RAND,sample) ;
    save_snps(name,snp_drivers,max_size,0,NULL) ; }
*/
/*
    delete [] snp_no ; delete [] snp_drivers ;
*/
/*
    if (nsam==1) {  // do this only when making images of tumours
      //printf("saving images...\n") ;
      int j=0 ;

      for (int i=0;i<genotypes.size();i++) {
        if (genotypes[i]!=NULL && genotypes[i]->number>0) genotypes[i]->index=j++ ;
      }

      vecd li1(1,-1,-0.3) ; // lighting direction
      float density ;

      cout<<"NUM="<<NUM<<"; max_size="<< max_size<<endl;
      sprintf(name,"%s/2d_image1_%d.dat",NUM,max_size) ; density=save_2d_image(name,li1) ;

      // if you want to save moreimages in a single run, add new lines like this
      //       vecd li2(1,-1,-1) ; sprintf(name,"%s/2d_image2_%d.dat",NUM,max_size) ; density=save_2d_image(name,li2) ;
      sprintf(name,"%s/cells_%d.dat",NUM,max_size) ; save_positions(name,1./density) ;

      sprintf(name,"%s/genotypes_%d.dat",NUM,max_size) ;SaveGenotypes(name) ;

      sprintf(name,"%s/most_abund_gens_%d.dat",NUM,max_size) ; save_most_abund_gens(name,most_abund) ;



       //if (TakeProbe(NUM,P)) num_adv+=1;




       //P.x_max=minx;
       //P.x_max=minx + 100;

       //sprintf(name,"%s/Probes.dat", NUM); TakeProbe(name,P,2);
    }
*/
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
    max_size    =  pt.get <int>
    ("Parameters.Evolution.stop_time_cell_num");
    driver_adv  =  pt.get <float>                                             //
    ("Parameters.Evolution.driver_adv");                                   //
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
  pars.probe_pars.del_area.del_x = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_x");
  pars.probe_pars.del_area.del_y = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_y");
  pars.probe_pars.del_area.del_z = pt.get<int>
    ("Parameters.ProbeComparisonParameters.Probe.DeltaArea.del_z");
  pars.etalon.piece_num = pt.get<int>
    ("Parameters.ProbeComparisonParameters.piece_num");
  pars.etalon.max_dist_probe_trees = pt.get<float>
    ("Parameters.ProbeComparisonParameters.max_dist_probe_trees");
  pars.etalon.threshold_frac_cell_diff = pt.get<float>
    ("Parameters.ProbeComparisonParameters.threshold_frac_cell_diff");
  pars.etalon.threshold_frac_mut_diff = pt.get<float>
    ("Parameters.ProbeComparisonParameters.threshold_frac_mut_diff");
  pars.threshold_simulation_num = pt.get <int>
    ("Parameters.threshold_simulation_num");

    growth0
    = 1;

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

boost::property_tree::ptree GetInitParams(char *name){


}

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
  if (s_setting_file_name.find("model") != std::string::npos) {
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

void CreateInferenceConfFile(Parameters pars) {
  string full_way_to_config_file = pars.catalog_file_name+"/"+
                                   "model_inference_info_check.info";

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
  pt.add("Parameters.ProbeComparisonParameters.piece_num",
          1);
  pt.add("Parameters.ProbeComparisonParameters.max_dist_probe_trees",
          1000);
  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_cell_diff",
          0.1);
  pt.add("Parameters.ProbeComparisonParameters.threshold_frac_mut_diff",
          0.2);

  pt.add("Parameters.threshold_simulation_num",
          pars.threshold_simulation_num);

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
          //CreateInferenceConfFile(pars);
//       // cout<<"analize frec()....."<<endl;
//       // AnalyzeFrec A(cells,L,100);

//        MainParameters main_pa/rs={driver_adv, gama, driver_prob, max_dist_probe_trees};
//
//        TakeProbe(NUM, ProbeArea, file_tree_comparison, main_pars);
//        //};
//        //cout<<"driver_adv="<< driver_adv;
//        //cout <<"num_adv="<<  num_adv<<"\n";
//        file_tree_comparison.close();


          //Initialize probe area (shoud be changed later)
//        SetEtalonProbe(pars);
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
//    init(ind);
//    //#pragma omp parallel for
//      //cout<<"Number of ex i="<<i<<endl;
//      //driver_adv=_drand48();
//      //cout<<"driver_adv="<<driver_adv<<endl;
//      num_adv=0;
//      //for(int m=0 ;m<10 ; m++) {
//      // cout<<"m="<<m<<endl;
//      init_min_max();
//      //cout<<"run()....."<<endl;
//      if (run(ind)){
//        PopulationRange pop_range;
//        init_PopulationRange(pop_range);
//        //P.x_max = maxx; P.x_min = maxx - 50; P.y_min=-5; P.z_min=-5;   P.y_max=5; P.z_max=5;
//        ofstream file_with_etalon_probe;
//        char name_file[256];
//        sprintf(name_file,"%s/file_tree_comparison.dat",NUM);
//        file_tree_comparison.open(name_file, ios::ate | ios::in);
//        if (file_tree_comparison.is_open()==false) file_tree_comparison.open(name_file);
//       // cout<<"TakeProbe()....."<<endl;
//        MainParameters main_pars={driver_adv, gama, driver_prob, max_dist_probe_trees};
//
//        TakeProbe(NUM, ProbeArea, file_tree_comparison, main_pars);
//        //};
//        //cout<<"driver_adv="<< driver_adv;
//        //cout <<"num_adv="<<  num_adv<<"\n";
//        file_tree_comparison.close();
//       // cout<<"analize frec()....."<<endl;
//
//       // AnalyzeFrec A(cells,L,100);
//      };
//
  };
  //end();
  //int i;
  //cin>>i;
  return 0;
}
