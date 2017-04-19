// Copyright 2017 Lebid Mykola all rights reserved


//Mykola begin
#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED
//Mykola end
#include <stdio.h>

#include <math.h>
#include <vector>
#include <string>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp



using namespace std;

double _drand48(void) ;
void _srand48(int a) ;
void err(char *reason) ;
void err(const char *reason) ;
void err(char *reason, int a) ;
void err(const char *reason, int a) ;
void err(char *reason, double a);
void err(const char *reason, double a);
void quicksort2(float *n, int *nums, int lower, int upper) ;
//Mykola begin
void init(bool is_for_etalon_sim, string string_catalog_name);
//Mykola end
int main_proc(int exit_size, int save_size, double max_time, double wait_time);
void end() ;
void reset() ;
void save_data() ;
void save_spatial(int *) ;
void save_snps(char *name,int *n, int total, int mode, int *) ;
float average_distance_ij() ;

extern int L ;
extern const int _resol, _bins ;
extern int volume, max_size ;

#ifndef classes_already_defined
#define classes_already_defined

#ifdef __linux
typedef unsigned int DWORD ;
typedef unsigned char BYTE ;
#else
#include <windows.h>
#endif

class vecd  // class of 3d vectors
{
public:
	double x, y, z;

	vecd(double x0, double y0, double z0) : x(x0), y(y0), z(z0) { }
	vecd() : x(0), y(0), z(0) {	}
	inline bool operator== (const vecd& V2) const
	 { return (x == V2.x && y == V2.y && z == V2.z) ; }
	inline vecd operator+ (const vecd& V2)
	 const { return vecd(x+V2.x,y+V2.y,z+V2.z) ; }
	inline vecd operator- (const vecd& V2)
	const { return vecd(x-V2.x,y-V2.y,z-V2.z) ; }
	inline vecd operator- () const { return vecd(-x,-y,-z);	}
	inline vecd operator/ (double S ) const {
		double f = 1.0/S;
		return vecd (x*f,y*f,z*f);
	}
	inline vecd operator* (const vecd& V2)
	const { return vecd (x*V2.x,y*V2.y,z*V2.z); }
	inline vecd operator* (double S) const { return vecd (x*S,y*S,z*S);	}
	inline void operator+= ( const vecd& V2 ) {
		x+=V2.x; y+=V2.y; z+=V2.z;
	}
	inline void operator-= ( const vecd& V2 ) {
		x -= V2.x; y -= V2.y; z -= V2.z;
	}
	inline void operator*= (double S) {
		x*=S ; y*=S ; z*=S ;
	}
	inline void operator/= (double S) {
    double S2=1./S ;
    x*=S2 ; y*=S2 ; z*=S2 ;
	}
};

inline double norm(vecd &a)
{
  return sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double squared(vecd &a)
{
  return (a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double scalar(vecd &a, vecd &b) {
  return a.x*b.x+a.y*b.y+a.z*b.z ;
}

inline vecd cross(vecd &a,vecd &b)
{
  vecd c ;
  c.x=a.y*b.z-a.z*b.y ;
  c.y=a.z*b.x-a.x*b.z ;
  c.z=a.x*b.y-a.y*b.x ;
  return c ;
}

inline void normalize(vecd &a)
{
  double l=sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
  a.x/=l ; a.y/=l ; a.z/=l ;
}

class IVec  // class of 3d integer vectors
{
  public:
	int i, j, k;
	IVec( int Ini, int Inj, int Ink ) : i( Ini ), j( Inj ), k( Ink ) { }
	IVec( ) : i(0), j(0), k(0) { }
	inline bool operator== (const IVec& V2) const
	{ return (i == V2.i && j == V2.j && k == V2.k); }
	inline void operator+= ( const IVec& V2 )
	{ i += V2.i; j += V2.j; k += V2.k; }
};


struct Cell {
  short unsigned int lesion ;
  short int x,y,z ;
  unsigned int gen;
};

#ifndef PUSHING
class Sites {
  public:
    DWORD *s ;
    Sites(int n0) { int n=1+(n0>>5) ; s=new DWORD[n] ;
    if (s==NULL) err("out of memory when allocating Sites") ;
    for (int i=0;i<n;i++) s[i]=0 ; }
    ~Sites() { delete [] s ; }
    inline void set(const unsigned int i) { s[(i>>5)]|=1<<(i&31) ; }
    inline void unset(const unsigned int i) { s[(i>>5)]&=~(1<<(i&31)) ; }
    inline int is_set(const unsigned int i) { return (s[(i>>5)]>>(i&31))&1 ; }
};
#else
typedef int Sites ;
#endif

extern vector <Cell> cells ;

#ifndef PUSHING
struct Lesion {
  int wx ;
  vecd r,rold,rinit ;
  double rad, rad0 ;
  int n,n0 ;
  vector <int> closest ;
  Sites **p ;
  static int nl ;
  static double maxdisp ;
  Lesion(int cellno, int g, int x0, int y0, int z0) {
    // cellno is not used anywhere in this version of the method
    rad=rad0=1 ;
    r=vecd(x0,y0,z0) ; rinit=rold=r ;
    closest.clear() ;
    wx=4 ; p=new Sites*[wx*wx] ;
    int i;
    for (i=0;i<wx*wx;i++) {
      p[i]=new Sites(wx) ;
    }
    Cell c ; c.x=c.y=c.z=0 ; c.gen=g ; c.lesion=nl++ ;
    if (nl>65000) err ("nl>65000") ;
    p[(wx/2)*wx+wx/2]->set(wx/2) ;
    cells.push_back(c) ; volume++ ; n=n0=1 ;
  }
  ~Lesion() {
    nl-- ;
    for (int i=0;i<wx*wx;i++) delete p[i] ;
    delete [] p ;
  }
  void update_wx() ;
  void find_closest() ;
  void one_move_step() ;
  void reduce_overlap() ;
  int no_free_sites(int x, int y, int z);
  void choose_nn(int &x, int &y, int &z);
};
#else
struct Lesion {
  int wx ;
  vecd r,rold,rinit ;
  double rad, rad0 ;
  int n,n0 ;
  vector <int> closest ;
  Sites ***p ;
  static int nl ;
  static double maxdisp ;
  Lesion(int cellno, int g, int x0, int y0, int z0) {
    rad=rad0=1 ;
    r=vecd(x0,y0,z0) ; rinit=rold=r ;
    closest.clear() ;
    wx=16 ; p=new Sites**[wx] ;
    int i,j,k;
    for (i=0;i<wx;i++) {
      p[i]=new Sites*[wx] ; if (p[i]==NULL) err("out of memory") ;
      for (j=0;j<wx;j++) {
        p[i][j]=new Sites[wx] ;
        for (k=0;k<wx;k++) p[i][j][k]=-1 ;
      }
    }
    Cell c ; c.x=c.y=c.z=0 ; c.gen=g ; c.lesion=nl++ ;
    if (nl>65000) err ("nl>65000") ;
    p[wx/2][wx/2][wx/2]=cellno ;
    cells.push_back(c) ; volume++ ; n=n0=1 ;
  }
  ~Lesion() {
    nl-- ;
    for (int i=0;i<wx;i++) {
      for (int j=0;j<wx;j++) delete p[i][j] ;
      delete [] p[i] ;
    }
    delete [] p ;
  }
  void update_wx() ;
  void find_closest() ;
  void one_move_step() ;
  void reduce_overlap() ;
  int no_free_sites(int x, int y, int z);
  // find direction of least drag
  void find_dir_min_drag(int i, int j, int k, int &in, int &jn, int &kn);
};
#endif

const unsigned int RESISTANT_PM = 1<<31 ;
const unsigned int DRIVER_PM = 1<<30 ;
const unsigned int L_PM = (1<<30) - 1 ;

struct Genotype {
  vector <unsigned int> sequence ;
  BYTE no_resistant, no_drivers ;
  // m = migration probability before/after treatment
  float death[2], growth[2], m[2] ;
  int number; // number of cells of this genotype total/on the surface
  int num_snp; // number of snp in genotype
  int prev_gen ;
  int index ; // this is set only when saving data
  Genotype(void) ;
  ~Genotype(void) { sequence.clear() ; }
  Genotype(Genotype *mother, int prevg, int no_snp) ;
};

class Hist {
  public:
  int x,n ;
  long long int x2 ;
  Hist() { x=n=0 ; x2=0 ; }
	inline void operator+= ( Hist& H2 ) { x+=H2.x ; x2+=H2.x2 ; n+=H2.n ; }
	inline void operator+= ( int val ) { x+=val ; x2+=val*val ; n++ ; }
  void r() { x=n=0 ; x2=0 ; }
};



#endif

extern vector<Genotype*> genotypes ;
extern vector<Lesion*> lesions ;

//Mykola begin
//Parameter descriptions
///////////////////////////////////////////////////////////////////////////////
struct Evolution {
    float driver_adv;
    float mutation_rate;
    float driver_mutation_rate;
    // boundary for number of cells in population
    int stop_time_cell_num;
};


struct DeltaArea {
    int del_x;
    int del_y;
    int del_z;
};
// structure ProbeArea
// Define Rectangular cuboid for probe cells
struct Area {
    int min_x;  int max_x;
    int min_y;  int max_y;
    int min_z;  int max_z;
};

// structure dynamic
struct Borders {
    int min_x;  int max_x;
    int min_y;  int max_y;
    int min_z;  int max_z;
};

struct Probe {
    bool is_random;
    Area area;
    DeltaArea del_area;
    Borders population_borders;
    int cell_num;
};

struct ComparisonResults {
    float dist = 10000;
    int num_mut = - 1;
    int e_num_mut = -1;
    int delta_mut = - 1;
};

struct ProbeComparisonParameters {
    string comparison_file_name;
    string e_mut_file_name;
    string cataloge_with_e_mut_file;
    int e_seed;
    Probe probe_pars;  // for etalon
    Evolution evo_pars;  // for etalon
    ComparisonResults comparison_results;
    unsigned int dist_type;
    float max_dist_probe_trees;
    float threshold_frac_cell_diff;  // threshold frac cell number discrepancy
    float threshold_frac_mut_diff;  // see above (for mutation)
    int piece_num; // number of probe pieces (attempts)
};


struct Parameters {
    string par_file_name;
    string file_name;
    string catalog_file_name;
    int seed;  // (RAND)
    int threshold_simulation_num;  // number of attempts to make simulation
                                   // that is important for small birth rate
    bool is_simulation_for_etalon_probe; // there is 0 - for trials
    Probe probe_pars;
    Evolution evo_pars;
    ProbeComparisonParameters etalon;
};
///////////////////////////////////////////////////////////////////////////////





//Mykola end


//Mykola begin
#endif // CLASSES_H_INCLUDED
//Mykola end
