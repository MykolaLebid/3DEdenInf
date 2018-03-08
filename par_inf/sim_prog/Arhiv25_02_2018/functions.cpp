#include <math.h>
#include "params.h"
#include "classes.h"
//#include "conf_structures.h"





int how_many_SNPs_identical(Genotype *a, Genotype *b)
{
  int i,j,n=0;
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++)
      if (a->sequence[i]==b->sequence[j]) n++ ;
  return n ;
}

int how_many_SNPs_identical(Genotype *a, Genotype *b, float cutoff, int *snp_no)
{
  int i,j,n=0, ntot=cells.size();
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++)
      if (a->sequence[i]==b->sequence[j] && snp_no[(a->sequence[i])&L_PM]>cutoff*ntot) n++ ;
  return n ;
}


void select_two_random_cells(int &i, int &j, int &n)
{
  int ntot=cells.size() ;
  double dist ;
  do {
    i=int(_drand48()*ntot) ; j=int(_drand48()*ntot) ;
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    dist=SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z) ;
  } while (_drand48()<0.5 && dist>_drand48()*100) ; // choose preferentially cells which are close to one another
  dist=sqrt(dist) ;

  n=int(dist/_resol) ;
  if (n<0 || n>=_bins) err("n",n);
}

float average_distance_ij()
{
  int ntot=cells.size() ;
  double avdist=0 ;
  for (int n=0;n<ntot;n++) {
    int i=int(_drand48()*ntot), j=int(_drand48()*ntot) ;
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    avdist+=sqrt(SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z)) ;
  }
  return (avdist/ntot) ;
}








