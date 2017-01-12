/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in

   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban,
   Bert Vogelstein, and Martin A. Nowak. “A Spatial Model Predicts That
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity.” Nature 525,
   no. 7568 (September 10, 2015): 261–64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/


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

#if (defined(GILLESPIE) && defined(FASTER_KMC)) || (defined(GILLESPIE) && defined(NORMAL)) || (defined(NORMAL) && defined(FASTER_KMC))
  #error too many methods defined!
#endif

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
extern int RAND, sample, treatment, max_size ;
extern double tt ;//time


float driver_adv;
float driver_prob;
float gama;

int num_adv;

int sample=0;



int minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;



void save_positions(char *name, float dz)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<1e4) fprintf(data,"%d %d %d %d %d\n",i,
                                                                     int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z),
                                                                     cells[i].gen);


        //fprintf(data,"%d %d %d %d %u\n",int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z), genotypes[cells[i].gen]->index) ;//WTF
  }
  fclose(data) ;
}


void search_min_max()
{
//Mykola search of min and max for all cells in all lesions
  int i;
  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  }
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


void save_genotypes(char *name)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;

    if (g!=NULL && g->number>0) {
      fprintf(data,"%d   %d   %d   (%d   %d)   %d\t"
              ,i,  g->prev_gen, g->num_snp, g->no_drivers, g->no_resistant, g->number) ;
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
      }
      fprintf(data,"\n");
    }
  }

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
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index) ;
    }
  }
  fclose(data) ;

}

bool run(int ind)
{
  int nsam=1;
  init(ind);
  //char name[256] ;
  //sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
  //FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) {
    reset();
#ifdef MAKE_TREATMENT_N
    int s=0 ; while (main_proc(max_size,-1,-1, 10)==1) { s++ ; reset() ; } ; // initial growth until max size is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
//    save_data() ;
    treatment=1 ;
    double max_time=2*tt ;
    main_proc(1.25*max_size,-1,max_time, 10) ; // treatment
#elif defined MAKE_TREATMENT_T
    cout<<"#elif defined MAKE_TREATMENT_T"<<endl;
    int s=0 ; while (main_proc(-1,-1,time_to_treat, 10)==1) { s++ ; reset() ; } ; // initial growth until max time is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
   // save_data() ;
    treatment=1 ;
    double max_time=2*tt ;
    max_size=cells.size()*1.25 ;
    main_proc(max_size,-1,max_time, 10) ; // treatment

#else
    int s=0 ; while (main_proc(max_size,2,-1, -1)==1 && (s<1000)) {
        //driver_adv=_drand48();
        s++ ; reset();
    } // initial growth until max size is reached

    //if (s>0) printf("resetted %d times\n",s);
    search_min_max();

    if (s<1000) return true;
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

      sprintf(name,"%s/genotypes_%d.dat",NUM,max_size) ; save_genotypes(name) ;

      sprintf(name,"%s/most_abund_gens_%d.dat",NUM,max_size) ; save_most_abund_gens(name,most_abund) ;



       //if (take_probe(NUM,P)) num_adv+=1;




       //P.x_max=minx;
       //P.x_max=minx + 100;

       //sprintf(name,"%s/Probes.dat",NUM); take_probe(name,P,2);
    }
*/
#endif
  }

}

void init_min_max()
{
    minx=1<<20, maxx=-minx, miny=minx, maxy=-minx, minz=minx, maxz=-minx;
}

void init_min_max_diff(Min_max_diff & A)
{
  A.min_x=minx;
  A.max_x=maxx;
  A.min_y=miny;
  A.max_y=maxy;
  A.min_z=minz;
  A.max_z=maxz;
  A.del_1=30;
  A.del_2=10;
  A.del_3=10;

}


int main(int argc, char *argv[])
{
#if defined(GILLESPIE)
  cout <<"method: GILLESPIE\n" ;
#endif
#if defined(FASTER_KMC)
  cout <<"method: FASTER_KMC\n" ;
#endif
#if defined(NORMAL)
  cout <<"method: NORMAL\n" ;
#endif



  //cout<<"save results of simulation as the etalon probe - '0'; compare with etalon probe - other symbol"<<endl;
  //int ind;cin>>ind;
  int ind=0;
  int nsam;
  float max_dist_for_prob;

  if (argc!=8) { err(" Error:: arguments needed: name, no_samples, RAND, driver_adv, driver_mutation_probability, mutation_probability, max dist for prob. Program terminated. \n"); }
  else {
    NUM=argv[1];
    nsam=atoi(argv[2]);
    RAND=atoi(argv[3]);
    driver_adv = atof(argv[4]);
    driver_prob= atof(argv[5]);
    gama= atof(argv[6]);
    max_dist_for_prob = atof(argv[7]);
  };
  //cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<" "<< driver_adv<<" "<<driver_prob<<" "<<gama<<" "<<max_dist_for_prob<<endl ;

  _srand48(RAND);


  Probe P;



  if (ind==1) {


    //Initialize main parameters
    driver_adv = 0.3;
    driver_prob = 0.00002;
    gama = 0.01;
    run(ind);

    //Initialize probe area (shoud be changed later)
    P.x_max = maxx; P.x_min = maxx - 30 ; P.y_min=-5; P.z_min=-5;  P.y_max=5; P.z_max=5;

    take_etalon_probe(NUM,P,20);


  } else {




        init(ind);
        //for (driver_adv=0.00; driver_adv<=1; driver_adv+=0.05) {
        //#pragma omp parallel for


        for (int i=0; i<1 ; i++){

            //cout<<"Number of ex i="<<i<<endl;
            //driver_adv=_drand48();
            //cout<<"driver_adv="<<driver_adv<<endl;
            num_adv=0;
            //for(int m=0 ;m<10 ; m++) {
            // cout<<"m="<<m<<endl;
            init_min_max();
            //cout<<"run()....."<<endl;
            if (run(ind)){

                Min_max_diff info;
                init_min_max_diff(info);

                //P.x_max = maxx; P.x_min = maxx - 50; P.y_min=-5; P.z_min=-5;   P.y_max=5; P.z_max=5;

                ofstream file_tree_comparison;
                char name_file[256];

                sprintf(name_file,"%s/file_tree_comparison.dat",NUM);
                file_tree_comparison.open(name_file, ios::ate | ios::in);

                if (file_tree_comparison.is_open()==false) file_tree_comparison.open(name_file);
               // cout<<"take_probe()....."<<endl;
                take_probe(NUM,info,file_tree_comparison, driver_adv, gama, driver_prob, max_dist_for_prob);

                //};
                //cout<<"driver_adv="<< driver_adv;
                //cout <<"num_adv="<<  num_adv<<"\n";
                file_tree_comparison.close();
               // cout<<"analize frec()....."<<endl;

               // Anal_frec A(cells,L,100);

            };

        };


  };

  //end();


  //int i;
  //cin>>i;

  return 0;
}
