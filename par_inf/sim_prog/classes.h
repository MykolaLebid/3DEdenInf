// Copyright 2017 Lebid Mykola all rights reserved
#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include <vector>
#include <string>
#include <deque>
#include <set>
#include <list>
#include "conf_structures.h" // structure Evolution for the function main_proc

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp



struct SimData;
struct SelectedGenotype;
class  Genotype;
struct Cell;
struct SelectedCell;

double _drand48(void) ;
void _srand48(int a) ;

void quicksort2(float *n, int *nums, int lower, int upper) ;
void init(bool is_for_etalon_sim, std::string string_catalog_name);

int mainProg3D(/*int exit_size, int save_size, double max_time,
							double wait_time,*/ Evolution & evolution, SimData & sim_data);

int mainProgWellMixed(Evolution & evolution, SimData & sim_data);

void cleanSimData(SimData & sim_data);

void initSimData(const float driver_adv, const float driver_mut_rate,
									SimData & sim_data);

#ifdef __linux
typedef unsigned int DWORD ;
typedef unsigned char BYTE ;
#else
//#include <windows.h>
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

//inline double scalar(vecd &a, vecd &b) {
//  return a.x*b.x+a.y*b.y+a.z*b.z ;
//}

//inline vecd cross(vecd &a,vecd &b)
//{
//  vecd c ;
//  c.x=a.y*b.z-a.z*b.y ;
//  c.y=a.z*b.x-a.x*b.z ;
//  c.z=a.x*b.y-a.y*b.x ;
//  return c ;
//}

//inline void normalize(vecd &a)
//{
//  double l=sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
//  a.x/=l ; a.y/=l ; a.z/=l ;
//}

//class IVec  // class of 3d integer vectors
//{
//  public:
//	int i, j, k;
//	IVec( int Ini, int Inj, int Ink ) : i( Ini ), j( Inj ), k( Ink ) { }
//	IVec( ) : i(0), j(0), k(0) { }
//	inline bool operator== (const IVec& V2) const
//	{ return (i == V2.i && j == V2.j && k == V2.k); }
//	inline void operator+= ( const IVec& V2 )
//	{ i += V2.i; j += V2.j; k += V2.k; }
//};


struct Cell{
//  short unsigned int lesion;
  short int x,y,z;
  Genotype * gen;
};

struct SelectedCell{
//short unsigned int lesion;
//short int x,y,z;
//	unsigned int gen;
Genotype * gen;
};

class Sites {
  public:
    unsigned int *s ;// 0 to 4,294,967,295
    Sites(int n0) {
    	unsigned int n = 1 + (n0>>5) ;
    	s = new unsigned int [n] ;
    // n_0 is a number of sites in 3D grid.
    // Ex. n_0 = 4 then we work with 4^3 sites in 1 array
    // n is dim of mass s (less then n in 32 times)
    // 1 - is occupied in s and 0 - is free
			if (s==NULL) {
				std::cout<<"out of memory when allocating Sites\n";
				exit(0);
			};
			for (unsigned int i = 0;	i < n;	i++) s[i] = 0;
    }

    ~Sites() {
    	delete [] s;
		}

    inline void set(const unsigned int i) { s[(i>>5)]|=1<<(i&31) ; } // mark i site as occupied
    inline void unset(const unsigned int i) { s[(i>>5)]&=~(1<<(i&31)) ; } // mark i site as free
    inline int is_set(const unsigned int i) { return (s[(i>>5)]>>(i&31))&1 ; } // is occupied?
};

//extern std::vector <Cell> cells ;

struct LesionCoords{// absolute coordinates of lesion location on a 3d grid
	int x;
	int y;
	int z;
};

class Lesion {// There is a centered dynamic 3D grid with border size border_size
public:
  unsigned int border_size = 4;
  vecd r, r_old, r_init;
  double rad, rad0;
  int n, n0;
//  static double maxdisp;
  std::vector <int> closest ;
  // the array of sites
  Sites **p ;
  // the total number of lesions in the moment of lesion appearance
  int current_num_lesions;
  Lesion(Cell& cell, int& total_num_lesions);
  ~Lesion();
  void update_border_size() ;
  void find_closest(std::vector<Lesion*> lesions) ;
  void one_move_step(std::vector<Lesion*> lesions) ;
//  void reduce_overlap(std::vector<Lesion*> lesions) ;
  int no_free_sites(int x, int y, int z);
  void choose_nn(int &x, int &y, int &z);
};

const unsigned int RESISTANT_PM = 1<<31 ;
const unsigned int DRIVER_PM = 1<<30 ;
const unsigned int L_PM = (1<<30) - 1 ;


struct SelectedGenotype{
	std::vector <unsigned int> sequence;
  unsigned int cell_num;					// number of cells of this genotype
	unsigned int gen_array_index;	//
	//	unsigned int current_num_snps;
};

struct GenotypeInfo{
  float birth;
  float death;
  unsigned int cell_number;  // number of cells of this genotype
};

class Genotype{
public:
  // Structural variables
  Genotype * ancestor_gen;
  std::set <Genotype *> descendant_gen_set;
	unsigned int gen_array_index; // for fast deletion of all genotypes

	//List of new SNPs in the genotype
  std::list <unsigned int> sequence;

  GenotypeInfo * genotype_info;
	//int no_drivers;

  //unsigned int num_snp; // number of snp in genotype
	//Genotype * index; 	// this is set only when saving data
	//unsigned int absolut_index;

/////////////////////////////////////////////////
/// \brief constructor of the first genotype
/////////////////////////////////////////////////
  Genotype(const float birth0, const float death0);
/////////////////////////////////////////////////
/// \brief construct of following genotypes
/// \param total_num_snps
///				 total number of single nucleotide polymorphisms
///				 on the moment of function call and
///				 it can be changed
/// \return
///
/////////////////////////////////////////////////
	Genotype(Genotype *	mother,
					 const unsigned int out_gen_array_index,
					 const float driver_adv,
					 const float driver_mutation_rate,
					 const unsigned int num_new_snps,
					 int & total_num_snps);
//support functions for genotype_info
	void incCellNum(){genotype_info->cell_number++;}
	void decCellNum(){genotype_info->cell_number--;}
	unsigned int getCellNum(){return genotype_info->cell_number;}
	float getBirth(){return genotype_info->birth;}
  float getDeath(){return genotype_info->death;}


  Genotype * getAncestorGen(){return ancestor_gen;}
  unsigned int getDescendantNum(){return descendant_gen_set.size();}
	void setAncestor(Genotype * gen){ancestor_gen = gen;}
	void addDescendant(Genotype * gen){descendant_gen_set.insert(gen);}


  void freeInfo();
  bool includeAncestorMuts();
	void delMotherConnection(){if(ancestor_gen!=nullptr)
																ancestor_gen->descendant_gen_set.erase(this);}
	~Genotype(void){
		freeInfo();
	};
};

struct SimData{
	std::vector<Genotype * > genotypes;
	std::vector<Lesion * > lesions;
	std::vector<Cell> cells;
	//Genotype * root_genotype = nullptr;
};


#endif // CLASSES_H_INCLUDED
