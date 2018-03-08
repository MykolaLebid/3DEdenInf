// Copyright 2017 Lebid Mykola all rights reserved
#include <cmath>

#include <vector>
#include <iostream>
#include <string>

#include "params.h"
#include "classes.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

//using namespace std;

//char *NUM ; // name given as 1st argument from the command line



#if defined __linux
#include <unistd.h>
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#elif defined __APPLE__
typedef unsigned int DWORD ;
int memory_taken()
{
  return 0 ; // not implemented
}
#else
#include <windows.h>
#include <psapi.h>
//Mykola begin
#include <time.h>
// Added to support GetProcessMemoryInfo()
//#pragma comment(lib, "psapi.lib")
//Mykola end
int memory_taken() // return memory available in MB
{
//	PROCESS_MEMORY_COUNTERS info;
//	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info));
//Mykola: add psapi with out -l to
//	return (int) (info.WorkingSetSize/(1<<20));
return 0;
}
//Mykola end

#endif

//void err(char *reason)
//{
//  std::cout <<reason<<endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//
//void err(const char *reason)
//{
//  std::cout <<reason<<endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}


template<typename Reason, typename Number>
void err(Reason s, Number num)
{
    std::cout << s << num << '\n';
#ifdef __WIN32
		system("pause") ;
#endif
		exit(0);
};

template<typename Reason>
void err(Reason s)
{
    std::cout << s << '\n';
#ifdef __WIN32
		system("pause") ;
#endif
		exit(0);
};


//
//void err(char *reason, int a) {
//	std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//
//void err(char *reason, unsigned int a) {
//	std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//void err(const char *reason, int a) {
//  std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//
//void err(char *reason, char *a)
//{
//  std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//void err(const char *reason, double a)
//{
//  std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}
//void err(char *reason, double a)
//{
//  std::cout <<reason<<": "<<a<<std::endl ;
//#ifdef __WIN32
//  system("pause") ;
//#endif
//  exit(0) ;
//}

static long long unsigned int _x=0x000100010001LL,
                              _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

void _srand48(int a) { _x=a ; }
//Mykola begin
void init(bool is_for_etalon_sim, std::string catalog_file_name);
//Mykola end
void end() ;

double /*tt=0,*/ tt_at_start ;
int start_clock ;

//int L=0 ; // total number of SNPs
int volume ; // total volume of the tumor
std::vector <int> drivers ; // vector of driver mutations
//FILE *drivers_file ;
int /*treatment=0,*/ cells_at_start;
FILE *times ;
extern int sample ;
//int RAND ; // random seed
char *timesbuffer ;

int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  const double l=exp(-gama) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;
}


std::vector<Cell> cells ;

//int Lesion::nl=0 ;
double Lesion::maxdisp=0 ;
double max_growth_rate ;

Genotype::Genotype(const float death0,
									 const float growth0) :
									 death(death0),
									 growth(growth0)
{
  number 			= 1;
  no_drivers	= 1;
  prev_gen		= - 1;

	sequence.clear();
};

Genotype::Genotype(	Genotype *mother,
									  const int mother_index,
										const unsigned int num_new_snps,
										unsigned int & total_num_snps) //
{
  death			 		 =	mother->death ;
  growth		 		 =	mother->growth;
  prev_gen	 		 =	mother_index;
  sequence	 		 =  mother->sequence;
  no_drivers 		 =  mother->no_drivers;


  for (unsigned int i = 0; i < num_new_snps; i++) {
    if ((driver_adv > 0 ) && (_drand48()<driver_prob)) {
			death *= 1 - driver_adv ;
      drivers.push_back(total_num_snps) ;
      sequence.push_back((total_num_snps++)|DRIVER_PM) ;
      no_drivers++;
    } else {
			sequence.push_back(total_num_snps++);
    };
  };


  if (total_num_snps>1e9) {
			std::cout<<"total number of SNPs is too big";
			exit(0);
	};

  number = 1;
}

std::vector<Genotype*> genotypes ;
std::vector<Lesion*> lesions ;

Lesion::Lesion(Cell& cell, int& total_num_lesions) {
		int x0 = cell.x;
		int y0 = cell.y;
		int z0 = cell.z;
    r  =  vecd(x0,y0,z0);

    r_init = r_old = r;
    closest.clear();
    border_size = 4;
    p=new Sites*[border_size*border_size];
    //int i;
    for (int i = 0; i<border_size*border_size; i++)
      p[i]=new Sites(border_size);
    //Cell c ;
    //c.x=c.y=c.z=0;
    //c.gen=g;
    //c.lesion=nl++;
    cell.x = 0;
    cell.y = 0;
    cell.z = 0;
    cell.lesion = total_num_lesions++;

    if (total_num_lesions>65000)
			err("total num lesions>65000");

    p[(border_size/2)*border_size+border_size/2]->set(border_size/2) ;
    //cells.push_back(c);
    volume++ ;
    n=n0=1 ;
};

Lesion::~Lesion()
{
	for (int i=0;i<border_size*border_size;i++) delete p[i] ;
	delete [] p ;
};


void Lesion::update_border_size()
{
  int i,j,k;
  int new_border_size=int(border_size*1.25) ;
  if (new_border_size%2==1) new_border_size++ ; // make sure it's even
  int diff=(new_border_size-border_size)/2 ;

#ifdef PUSHING
// fill in -1 in new sites  and extant lattice
  for (i=0;i<border_size;i++) {
    for (j=0;j<border_size;j++) {
      Sites *np=new Sites[new_border_size] ;
      for (k=0;k<new_border_size;k++) np[k]=-1 ;
      for (k=diff;k<border_size+diff;k++) np[k]=p[i][j][k-diff] ;
      delete p[i][j] ;
      p[i][j]=np ;
    }
  }

  Sites ***np=new Sites**[new_border_size] ;
  for (i=0;i<new_border_size;i++) np[i]=new Sites*[new_border_size] ;

  for (i=0;i<new_border_size;i++) {
    for (j=0;j<new_border_size;j++) {
      if (i<diff || i>=border_size+diff || j<diff || j>=border_size+diff) {
        np[i][j]=new Sites[new_border_size] ;
        for (k=0;k<new_border_size;k++) np[i][j][k]=-1 ;
      } else {
        np[i][j]=p[i-diff][j-diff] ;
      }
    }
  }

  for (i=0;i<border_size;i++) delete [] p[i] ;

#else
// extant lattice
  for (i=0;i<border_size*border_size;i++) {
    Sites *np=new Sites(new_border_size) ;
    for (k=diff;k<border_size+diff;k++) if (p[i]->is_set(k-diff)) np->set(k) ;
    delete p[i] ;
    p[i]=np ;
  }

  if (float(new_border_size)*float(new_border_size)>2e9)
		err("new border size too large",new_border_size) ;
  Sites **np=new Sites*[new_border_size*new_border_size] ;

  for (i=0;i<new_border_size;i++) {
    for (j=0;j<new_border_size;j++) {
      if (i<diff || i>=border_size+diff || j<diff || j>=border_size+diff) {
        np[i*new_border_size+j]=new Sites(new_border_size) ;
      } else {
        np[i*new_border_size+j]=p[(i-diff)*border_size+j-diff] ;
      }
    }
  }
#endif

  delete [] p ;
  p=np ;
  border_size=new_border_size ;

}

void Lesion::one_move_step() {
  /*int i,j;*/
  double mthis=this->n ;
  for (unsigned int i=0;i<closest.size();i++) {
    vecd dr=lesions[closest[i]]->r - this->r ;
    double r2=squared(dr), sumrad2=SQR(this->rad+lesions[closest[i]]->rad) ;
    if (r2<sumrad2) {
      double mi=lesions[closest[i]]->n ;
      double disp=(sqrt(sumrad2/r2)-1) ;
      if (fabs(disp)>maxdisp) maxdisp=fabs(disp) ;
      dr*=disp*1.1 ;
      this->r-=dr*mi/(mi+mthis) ;
      lesions[closest[i]]->r+=dr*mthis/(mi+mthis) ;
    }
  }
}

void Lesion::find_closest()
{
  r_old = r;
  closest.clear() ;
  for (unsigned int i=0; i<lesions.size(); i++) {
    vecd dr=this->r - lesions[i]->r ;
    double r2=squared(dr) ;
    if (r2>0 && r2<2*(SQR(this->rad+lesions[i]->rad))) {
      closest.push_back(i) ;
    }
  }
}

void Lesion::reduce_overlap()
{
  int i,k,temp;
  int *ind=new int[lesions.size()] ;
  for (unsigned int j=0;j<lesions.size();j++) ind[j]=j ;
  do {
    maxdisp=0 ;
    for (unsigned int j=0;j<lesions.size();j++)
			{ k=_drand48()*lesions.size() ; SWAP(ind[j],ind[k]) ; }
    for (unsigned int j=0;j<lesions.size();j++) {  // go through a random permutation
      i=ind[j] ;
      lesions[i]->one_move_step() ;

      vecd dr=lesions[i]->r - lesions[i]->r_old ;
      if (squared(dr)>SQR(lesions[i]->rad)) lesions[i]->find_closest() ;
    }
  } while (maxdisp>1e-2) ;
  delete [] ind ;
}

void reset(const int driver_adv, const int driver_mut_rate)
{
  /*tt=0 ; L=0 ;*/

//genotypes preparation
	for (auto &i:genotypes)
		if (i!=NULL) delete i;

//  for (unsigned int i=0;i<genotypes.size();i++)
//		if (genotypes[i]!=NULL) delete genotypes[i];
  genotypes.clear();
	const float growth0 = 1;
	const float death0 = 1 - driver_adv;
  max_growth_rate = growth0;
  genotypes.push_back(new Genotype(growth0,death0));

// lesion preparation
	for (auto &i:lesions)
		delete i;
//  for (unsigned int i=0; i<lesions.size();i++)
//		delete lesions[i];
	lesions.clear();
  cells.clear();
  drivers.clear();

  volume = 0;
  Cell cell;
  cell.x=cell.y=cell.z=0;
  cell.gen=0;

  int total_num_lesions = 0;

  lesions.push_back(new Lesion(cell, total_num_lesions));
  cells.push_back(cell);
}

#ifdef MOORE_NEIGHBOURHOOD
const int _nonn=26 ;
const int kx[27]={0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1},
          ky[27]={0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1},
          kz[27]={0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1};
int kln[27] ; // this is filled with lengths of (kx,ky,kz)
#endif

#ifdef VON_NEUMANN_NEIGHBOURHOOD
const int _nonn=6 ;
const int kx[7]={0,1,-1,0,0,0,0},
          ky[7]={0,0,0,1,-1,0,0},
          kz[7]={0,0,0,0,0,1,-1};
int kln[7] ; // this is filled with lengths of (kx,ky,kz)
#endif
//Mykola begin
//#if defined __linux
//void init(bool is_for_etalon_sim, string catalog_name) {
//  int i;
//  for (i=0;i<=_nonn;i++) kln[i]=sqrt(1.*SQR(kx[i])+1.*SQR(ky[i])+1.*SQR(kz[i])) ;
//Mykola begin
//  if (is_for_etalon_sim){
//    char txt[256] ;
//    sprintf(txt,"mkdir %s",catalog_name.c_str()) ; system(txt);
//  };
//  start_clock=clock() ;
//}
//#elif defined __APPLE__
//  // not defined yet
//#else
void init(bool is_for_etalon_sim, std::string path)
{
  int i;
  for (i=0;i<=_nonn;i++)
    kln[i]=sqrt(1.*SQR(kx[i])+1.*SQR(ky[i])+1.*SQR(kz[i])) ;
  if (is_for_etalon_sim){
    char txt[256];
    sprintf(txt,"mkdir %s",path.c_str()); system(txt);
  };
  start_clock = clock();
}
//#endif
//Mykola end

//void end() {
//  fclose(times) ;
//}

#ifdef PUSHING
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++) nfree-=p[(border_size+z+kz[n])%border_size]
		 [(border_size+y+ky[n])%border_size][(border_size+x+kx[n])%border_size]==-1?0:1 ;
  return nfree ;
}
#else
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++)
    nfree-=p[((border_size+z+kz[n])%border_size)*border_size +
					 (border_size+y+ky[n])%border_size]->
																	is_set((border_size+x+kx[n])%border_size) ;
  return nfree ;
}
inline void Lesion::choose_nn(int &x, int &y, int &z)
{
  static int nns[_nonn] ;
  int no=0,n ;
  for (n=1;n<=_nonn;n++)
    if (p[((border_size+z+kz[n])%border_size)*border_size +
				(border_size+y+ky[n])%border_size]->
				   is_set((border_size+x+kx[n])%border_size)==0) nns[no++]=n;
  if (no==0) { x=-1000000 ; return ; }
  n=nns[int(_drand48()*no)] ;
  z=(border_size+z+kz[n])%border_size ;
  y=(border_size+y+ky[n])%border_size ; x=(border_size+x+kx[n])%border_size;
}
#endif


inline int free_sites(int n)
{
  Lesion *ll=lesions[cells[n].lesion] ;
  int border_size=ll->border_size ;
  int k=cells[n].x+border_size/2, j=cells[n].y+border_size/2, i=cells[n].z+border_size/2 ;
  return ll->no_free_sites(k,j,i) ;
}


void quicksort2(float *n, int *nums, int lower, int upper)
{
	int i, m, temp ;
  float pivot, tempd;

	if (lower < upper)
	{
		SWAPD(n[lower], n[(upper + lower) / 2]); SWAP(nums[lower], nums[(upper + lower) / 2]);
		pivot = n[lower];
		m = lower;
		for (i = lower + 1; i <= upper; i++)
			if (n[i] > pivot)
			{
				m++;
				SWAPD(n[m], n[i]); SWAP(nums[m], nums[i]);
			}
		SWAPD(n[lower], n[m]); SWAP(nums[lower], nums[m]);
		quicksort2(n, nums, lower, m - 1);
		quicksort2(n, nums, m + 1, upper);
	}
}



#ifdef PUSHING
void Lesion::find_dir_min_drag(int i, int j, int k, int &in, int &jn, int &kn)       // find direction of least drag
{
  int nn, in0,jn0,kn0 ;
  std::vector <IVec> vis ; // vector of visited sites. it will begin with (i,j,k) and end at empty site

rep:
  vis.clear() ;
  vis.push_back(IVec(i,j,k)) ;
  in0=i ; jn0=j ; kn0=k ;

  do {      // loop goes over subsequent pushing events
    float mind=wx ;
    nn=-1 ;
    for (int nnnn=0;nnnn<10;nnnn++) {
      int nnn=1+_drand48()*_nonn ;
      in=in0 ; jn=jn0 ; kn=kn0 ;
      for (float drag=0;drag<mind;drag+=kln[nnn]) {
        in+=kz[nnn] ; jn+=ky[nnn] ; kn+=kx[nnn] ;
        for (int vv=0;vv<vis.size();vv++) if (vis[vv]==IVec(in,jn,kn)) goto brk ; // reject if trajectory passes through prev. visited sites
        if (p[in][jn][kn]==-1) { mind=drag ; nn=nnn ; break ; }
      }
brk:  continue ;
    }
    if (nn==-1) goto rep ;
        // now nn gives the direction of pushing

    in0+=kz[nn] ; jn0+=ky[nn] ; kn0+=kx[nn] ; // update position of the cell to be pushed
    vis.push_back(IVec(in0,jn0,kn0)) ; // and remember it...
  } while (p[in0][jn0][kn0]!=-1) ; // if the next position contains an empty site then exit

  // push all remembered cells except mother to make space for a single new daughter cell
  Sites sup=-1 ;        // sup is the elevated cell that needs to be inserted into new position
  for (int vv=1;vv<vis.size();vv++) {
    in0=vis[vv].i ; jn0=vis[vv].j ; kn0=vis[vv].k ;
    Sites snew=p[in0][jn0][kn0] ;
    p[in0][jn0][kn0]=sup ;
    if (sup!=-1) {
      Cell *c=&cells[sup] ;
      c->x=kn0-wx/2 ; c->y=jn0-wx/2 ; c->z=in0-wx/2 ;
    }
    sup=snew ;
  }

  in=vis[1].i ; jn=vis[1].j ; kn=vis[1].k ; // the new cell will be the second position from the list (1st is the mother cell which is not pushed)
  if (p[in][jn][kn]!=-1) err("!!!") ;
}
#endif


//-----------------------------------------------------
#if defined(NORMAL)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  unsigned int n; // index of a randomly chosen cell
  double time = 0; // timer of the evolution
  int /*i,j,k,*/ntot;
  short int total_number_lesions = 1;
  unsigned int genotype_death = 1;
  // Total number of single nucleotide polymorphisms.
  // Can be increased during a sell devision (creation)
  // within new Genotype initiation: genotypes.push_back
  unsigned int total_num_snps = 0;
////////////////////////////////////////////////////////////////////////////////
////description
////////////////////////////////////////////////////////////////////////////////
////  main_loop{
////  	time is change (according to tsc times number of cells ) and balance in step of birth event
////    take a cell randomly (n)
////    birth event is taken with probability witch is balanced with respect to time change tsc*genotypes[cells[n].gen]->growth[treatment]
////    {
////
////			if defined(CONST_BIRTH_RATE) take one free neighbor (equal probabilities)
////      if not defined(CONST_BIRTH_RATE) take one neighbor (equal probabilities). If it is free then  create cell and if not then don't create
////      if defined(VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC) take try the second time and take one free neighbor (equal probabilities)
////		  new and old cells will get mutations proportionally to two Poisson distributions with the same intensity
////      after then old cell takes death event with probability witch is balanced with respect to time change
////      tsc*genotypes[cells[n].gen]->death[treatment]
////    };
////
////
////
////
////  };
////////////////////////////////////////////////////////////////////////////////
  for(;;) {      // main loop
    //double tsc = 0.01*cells.size() ; if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ;
    time+=/*tsc**/timescale/cells.size() ;
    n=_drand48()*cells.size() ;
    Lesion *ll=lesions[0/*cells[n].lesion*/] ;
    //size of grid (wx - from -wx/2 to wx/2)
    unsigned int border_size = ll->border_size;
    int k=cells[n].x+border_size/2 ;
    int j=cells[n].y+border_size/2 ;
    int i=cells[n].z+border_size/2 ;
    int need_border_size_update=0 ;
		//if grid is too small change it
    if (k<2 || k>=((int)border_size - 3) ||
				j<2 || j>=((int)border_size - 3) ||
				i<2 || i>=((int)border_size - 3)) need_border_size_update=1 ;
#ifdef PUSHING
    if (ll->p[i][j][k]!=n) err("ll->p[i][j][k]!=n, p=",ll->p[i][j][k]) ;
#else
    if (ll->p[i*border_size+j]->is_set(k)==0) err("ll->p[i][j][k]==0, border size=",border_size) ;
#endif
    // cell reproduces if (tsc or 1/max_growth_rate ) * growth rate  (rebalance)
    // (birth rate don't grow up for standard case)
    if (_drand48()</*tsc**/genotypes[cells[n].gen]->growth) { // reproduction (treatment = 0) and growth will be as "normal"
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING) //begin***
		int nn = 1 + int(_drand48()*_nonn);//choose 1 neighbor
		int in=(border_size+i+kz[nn])%border_size,
		    jn=(border_size+j+ky[nn])%border_size,
		    kn=(border_size+k+kx[nn])%border_size;
#ifdef VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC // one more trial to find an empty site
		if (ll->p[in*border_size+jn]->is_set(kn)==1) {
        nn=1+int(_drand48()*_nonn) ;
        in=(border_size+i+kz[nn])%border_size ;
        jn=(border_size+j+ky[nn])%border_size ;
        kn=(border_size+k+kx[nn])%border_size ;
      }
#endif // VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC
      if (ll->p[in*border_size+jn]->is_set(kn)==0) {
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
				int in=i, jn=j, kn=k ;
				ll->choose_nn(kn,jn,in) ;
				if (kn!=-1000000) { // if there is at least one empty n.n., then.....
#elif defined(PUSHING)
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;
      {
#else
  #error inconsistent growth conditions
#endif //end***
        int no_SNPs=poisson() ; // newly produced cell mutants
//        if (_drand48()>genotypes[cells[n].gen]->m[treatment]) {
          Cell c ;
          c.x=kn-border_size/2 ;
          c.y=jn-border_size/2 ;
          c.z=in-border_size/2 ;
          c.lesion=cells[n].lesion ;
#ifndef PUSHING
          ll->p[in*border_size+jn]->set(kn) ; // set grid occupied
#else
          ll->p[in][jn][kn]=cells.size() ; // set grid occupied with number of cell
#endif
          if (no_SNPs>0) {// if number of SNP is not 0, then create a new genotype
            c.gen=genotypes.size() ;
            genotypes.push_back(new Genotype(genotypes[cells[n].gen],
																						 cells[n].gen,
																						 no_SNPs,
																						 total_num_snps)); // add new genotype
          } else { // if number of SNP is 0 add 1 cell in a correspondent genotype
            c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ;
          }
          cells.push_back(c) ; volume++ ; // add the new cell
          ll->n++ ;
#ifndef NO_MECHANICS
          //collision of two balls
          double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
          if (ll->rad/ll->rad0>1.05) {
            ll->reduce_overlap() ;
            ll->find_closest() ;
            ll->rad0=ll->rad ;
            ll->n0=ll->n ;
          }
#endif
//        } else { // make a new lesion (in our case 0 chance)
//          int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
//          if (no_SNPs>0) {
//            genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
//            lesions.push_back(new Lesion(cells.size(),genotypes.size()-1,x,y,z)) ;
//          } else {
//            genotypes[cells[n].gen]->number++ ;
//            lesions.push_back(new Lesion(cells.size(),cells[n].gen,x,y,z)) ;
//          }
//#ifndef NO_MECHANICS
//          lesions[lesions.size()-1]->find_closest() ;
//#endif
//        }
// BOTH_MUTATE
        no_SNPs=poisson() ; // old cell mutates
        if (no_SNPs>0) {
          genotypes[cells[n].gen]->number-- ;
          /*int pn=genotypes.size() ;*/
          genotypes.push_back(new Genotype(genotypes[cells[n].gen],
																					 cells[n].gen,
																					 no_SNPs,
																					 total_num_snps)) ;
          cells[n].gen=genotypes.size()-1 ;
          if (genotypes[cells[n].gen]->number<=0) {
            delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
          }
        }
      }
    }
#ifdef CORE_IS_DEAD //delete core cells but leave a cite occupied
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) {
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ;
    }
#endif

// now we implement death
//    if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ; // this is an alternative way but it does not change anything so no need to use it
#ifdef DEATH_ON_SURFACE
    if (genotypes[cells[n].gen]->death>0 && _drand48()<tsc*genotypes[cells[n].gen]->death*ll->no_free_sites(k,j,i)/float(_nonn))  { // death on the surface
#else
    if (_drand48()</*tsc**/genotypes[cells[n].gen]->death) { // death in volume
#endif
#ifndef PUSHING
      ll->p[i*border_size+j]->unset(k) ;
#else
      ll->p[i][j][k]=-1 ;
#endif
      ll->n-- ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ;
        for (unsigned int i=0;i<border_size;i++)
					for (unsigned int j=0;j<border_size;j++)
						for (unsigned int k=0;k<border_size;k++) {
#ifndef PUSHING
          double d=SQR(i-border_size/2)+SQR(j-border_size/2)+SQR(k-border_size/2) ; if (ll->p[i*border_size+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ;
#else
          double d=SQR(i-border_size/2)+SQR(j-border_size/2)+SQR(k-border_size/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ;
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        unsigned int nn=cells[n].lesion ;
        ll=NULL;
        delete lesions[nn] ;
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (unsigned int i=0;i<cells.size();i++)
						if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ;
        }
        lesions.pop_back();
        total_number_lesions--;
#ifndef NO_MECHANICS
        for (unsigned int i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;
        }
#endif
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) {
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
        genotype_death++;
      }
      if (n!=cells.size()-1) {
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->border_size/2, jj=cells[n].y+ll2->border_size/2, kk=cells[n].x+ll2->border_size/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_border_size_update && ll!=NULL) ll->update_border_size() ;

#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif
//Mykola begin
  //  if (wait_time>0 && time>tt_old+wait_time) { tt_old=time ; save_data(); }
  //  if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }
//Mykola end
    if (cells.size()==0) return 1 ;
    if (max_time>0 && time>max_time) return 3 ;
    if (exit_size>0 && ntot>=exit_size) {
				std::cout<< ntot<<"-"<< genotypes.size() <<" - "<<genotype_death<<"-"<<total_num_snps<<"\n";
				return 4 ;
		};

  }

}
#endif // NORMAL


//-----------------------------------------------------------------------------

#if defined(FASTER_KMC) || defined(GILLESPIE)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j,k,n,l,in,jn,kn,ntot;
  int cc=0, timeout=0 ;
	double tt = 0;

  for(;;) {      // main loop

#ifdef GILLESPIE //beginGILLESPIE
  // Gillespie
    double tot_rate = 0;
    //calculate total rate
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[0] * free_sites(n)/float(_nonn) ; // according to B Eden model
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) tot_rate+=genotypes[cells[n].gen]->growth[0] ; // according to C Eden model
#elif defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[0] ; // according to C article pushing case
#else
      #error inconsistent growth conditions
#endif
#ifdef DEATH_ON_SURFACE
      tot_rate+=genotypes[cells[n].gen]->death[0] * free_sites(n)/float(_nonn) ; // death on the surface
#else
      tot_rate+=genotypes[cells[n].gen]->death[0] ;  // death in volume
#endif
    }
    tt+=-log(1-_drand48())*timescale/tot_rate ;
    double q=_drand48()*tot_rate,r=0 ;
    int mode=0 ;
    //choose an event: birth -> mode = 0; death -> mode = 1
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[0] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) r+=genotypes[cells[n].gen]->growth[0] ;
#elif defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[0] ;
#else
      #error inconsistent growth conditions
#endif
      if (r>q) break ;
#ifdef DEATH_ON_SURFACE
      r+=genotypes[cells[n].gen]->death[0] * free_sites(n)/float(_nonn)  ; // death on the surface
#else
      r+=genotypes[cells[n].gen]->death[0] ;  // death in volume
#endif
      if (r>q) { mode=1 ; break ; }
    }
    if (n==cells.size()) err("n==cells.size() at t=",tt) ;
#endif //endGILLESPIE

#ifdef FASTER_KMC //beginFASTER_KMC
    double max_death_rate=1 ;
    double tot_rate=cells.size()*(max_growth_rate+max_death_rate) ;
    tt+=-log(1-_drand48())*timescale/tot_rate ;
    n=_drand48()*cells.size() ;
    double q=_drand48()*(max_growth_rate+max_death_rate), br,dr  ;
    int mode=0 ;
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING) //begin***
    br=genotypes[cells[n].gen]->growth[0] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
    if (free_sites(n)>0) br=genotypes[cells[n].gen]->growth[0] ; else br=0 ;
#elif defined(PUSHING)
    br=genotypes[cells[n].gen]->growth[0] ;
#else
      #error inconsistent growth conditions
#endif //end***
#ifdef DEATH_ON_SURFACE //begin****
    dr=genotypes[cells[n].gen]->death[0] * free_sites(n)/float(_nonn) ; // death on the surface
#else
    dr=genotypes[cells[n].gen]->death[0] ;  // death in volume
#endif//end****
    if (q<br) mode=0 ;
    else if (q<br+dr) mode=1 ;
    else mode=2 ;
#endif //endFASTER_KMC

    Lesion *ll=lesions[cells[n].lesion] ;
    int wx=ll->wx ;
    k=cells[n].x+wx/2 ; j=cells[n].y+wx/2 ; i=cells[n].z+wx/2 ;
    int need_wx_update=0 ;
    if (k<2 || k>=ll->wx-3 || j<2 || j>=ll->wx-3 || i<2 || i>=ll->wx-3) need_wx_update=1 ;

    if (mode==0) { // reproduction
#if !defined(PUSHING)
      int in=i, jn=j, kn=k ;
      ll->choose_nn(kn,jn,in) ;
#else
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;
#endif

      int no_SNPs=poisson() ; // newly produced cell mutants

			Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
#ifdef PUSHING
        ll->p[in][jn][kn]=cells.size() ;
#else
        ll->p[in*wx+jn]->set(kn) ;
#endif
        if (no_SNPs>0) {
          c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate
        } else {
          c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ;
        }
        cells.push_back(c) ; volume++ ;

        ll->n++ ;
#ifndef NO_MECHANICS
        double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
        if (ll->rad/ll->rad0>1.05) {
          ll->reduce_overlap() ;
          ll->find_closest() ;
          ll->rad0=ll->rad ;
          ll->n0=ll->n ;
        }
#endif

// BOTH_MUTATE
      no_SNPs=poisson() ; // old cell mutates
      if (no_SNPs>0) {
        genotypes[cells[n].gen]->number-- ;
        int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
        cells[n].gen=genotypes.size()-1 ;
        if (genotypes[cells[n].gen]->number<=0) {
          delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
        }
      }
    } //end_reproduction

#ifdef CORE_IS_DEAD
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) {
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ;
    }
#endif

// now we implement death
    if (mode==1) {
#ifdef PUSHING
      ll->p[i][j][k]=-1 ;
#else
      ll->p[i*wx+j]->unset(k) ;
#endif
      ll->n-- ;
//      if (ll->n<0) err("ll->n<0") ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ;
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
#ifdef PUSHING
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ;
#else
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ;
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        int nn=cells[n].lesion ;
//        if (ll!=lesions[nn]) err("ll!") ;
        ll=NULL ;
        delete lesions[nn] ;
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ;
        }
        lesions.pop_back() ;
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;
        }
#endif
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) {
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ;
      }
      if (n!=cells.size()-1) {
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->wx/2, jj=cells[n].y+ll2->wx/2, kk=cells[n].x+ll2->wx/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_wx_update && ll!=NULL) ll->update_border_size();

#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif
//Mykola begin
  //  if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
  //  if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }
//Mykola end
    if (cells.size()==0) return 1 ;
    if (max_time>0 && tt>max_time) return 3 ;
    if (exit_size>0 && ntot>=exit_size) return 4 ;

  }

}

#endif // FASTER_KMC or GILLESPIE
