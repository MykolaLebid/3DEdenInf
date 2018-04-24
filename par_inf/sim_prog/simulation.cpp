// Copyright 2017 Lebid Mykola all rights reserved
#include <cmath>

#include <vector>
#include <iostream>
#include <string>
#include <set>

#include "params.h"
#include "classes.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

//std::vector<Genotype*> genotypes ;
//std::vector<Lesion*> lesions ;



//using namespace std;

//char *NUM ; // name given as 1st argument from the command line

//#if defined __linux
//#include <unistd.h>
//typedef unsigned int DWORD ;
//#elif defined __APPLE__
//typedef unsigned int DWORD ;
//#else
//#endif
//
//#if defined __linux
//#include <unistd.h>
//typedef unsigned int DWORD ;
////int memory_taken() // return memory available in MB
////{
////  long long int rss = 0L;
////	FILE* fp = NULL;
////	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
////		return (size_t)0L;		/* Can't open? */
////	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
////	{
////		fclose( fp );
////		return (size_t)0L;		/* Can't read? */
////	}
////	fclose( fp );
////	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
////	return (int) (mem/(1<<20));
////}
////#include <sys/sysinfo.h>
////unsigned int freemem() // returns available memory in MB
////{
////  struct sysinfo s ;
////  sysinfo(&s) ;
////  return ((s.freeram)>>20) ;
////}
//#elif defined __APPLE__
//typedef unsigned int DWORD ;
////int memory_taken()
////{
////  return 0 ; // not implemented
////}
//#else
////#include <windows.h>
////#include <psapi.h>
////Mykola begin
////#include <time.h>
//// Added to support GetProcessMemoryInfo()
////#pragma comment(lib, "psapi.lib")
////Mykola end
////int memory_taken() // return memory available in MB
////{
//////	PROCESS_MEMORY_COUNTERS info;
//////	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info));
//////Mykola: add psapi with out -l to
//////	return (int) (info.WorkingSetSize/(1<<20));
////return 0;
////}
////Mykola end
//
//#endif
//
////void err(char *reason)
////{
////  std::cout <<reason<<endl ;
////#ifdef __WIN32
////  system("pause") ;
////#endif
////  exit(0) ;
////}
////
////void err(const char *reason)
////{
////  std::cout <<reason<<endl ;
////#ifdef __WIN32
////  system("pause") ;
////#endif
////  exit(0) ;
////}


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
void init(bool is_for_etalon_sim, std::string catalog_file_name);
void end() ;

//double /*tt=0,*/ tt_at_start ;
//int start_clock ;

//int L=0 ; // total number of SNPs
//int volume ; // total volume of the tumor
//std::vector <int> drivers ; // vector of driver mutations
//FILE *drivers_file ;
//int /*treatment=0,*/ cells_at_start;
//FILE *times ;
//extern int sample ;
//int RAND ; // random seed
//char *timesbuffer ;

int poisson(const float mut_rate)  // generates k from P(k)=exp(-mut_rate) mut_rate^k / k!
{
  const double l=exp(-mut_rate) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;
}


//std::vector<Cell> cells ;
//int Lesion::nl=0 ;
//double Lesion::maxdisp=0 ;
//double max_growth_rate ;

Genotype::Genotype(/*const float birth0,*/
									 const float death0)
{
	absolute_index = 0;
  //birth = 1;
  death = death0;
	cell_number = 1;
	descendant_number = 0;
  snps_number = 1;
  gen_array_index	= 0;
  ancestor_gen		= nullptr;
};

Genotype::Genotype(	Genotype *					mother,
										const unsigned int	out_gen_array_index,
										const float 				driver_adv,
										const float 				driver_mutation_rate,
										const unsigned int	num_new_snps,
										int &			total_num_snps)
{
  ancestor_gen	= mother;
	mother->descendant_number++;

	absolute_index	=	total_num_snps + 1;
	descendant_number = 0;

  //birth = mother->birth;
  death = mother->death;
	cell_number = 1;
  gen_array_index	= out_gen_array_index;
	snps_number = 0;
  for (unsigned int i = 0; i < num_new_snps; i++) {
    if ((driver_adv > 0 ) && (_drand48()<driver_mutation_rate)) {
			death *= 1 - driver_adv;
//      drivers.push_back(total_num_snps) ;
//      sequence.push_back((total_num_snps++)|DRIVER_PM) ;
//      no_drivers++;
    }
//     else {
//			sequence.push_back(total_num_snps++);
//    };
		snps_number++;
  };
	total_num_snps += snps_number;
  if (total_num_snps>1e9) {
			std::cout<<"total number of SNPs is too big";
			exit(0);
	};
}

//void Genotype::freeInfo()
//{
//	if (genotype_info != nullptr){
//		delete genotype_info;
//		genotype_info = nullptr;
//	};
////	for(auto &i: descendant_gen_vec){
////		i->ancestor_gen = this->ancestor_gen;
////		i->sequence.insert(i->sequence.begin(),
////											 this->sequence.begin(),
////											 this->sequence.end());
////	};
//};
bool Genotype::includeAncestorMuts()
{
	// cut down ancestor
	//check that ancestor_gen!=nullptr
	if ((ancestor_gen->getDescendantNum() == 1) && (ancestor_gen->getCellNum()==0)) {
//		sequence.splice(sequence.begin(), ancestor_gen->sequence);
		snps_number += ancestor_gen->snps_number;
		return true;
	} else {
		return false;
	};
};

//void Genotype::delMotherConnection()
//{
//};


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
    for (unsigned int i = 0; i < border_size * border_size; i++)
      p[i]=new Sites(border_size);
    //Cell c ;
    //c.x=c.y=c.z=0;
    //c.gen=g;
    //c.lesion=nl++;
    cell.x = 0;
    cell.y = 0;
    cell.z = 0;
    //cell.lesion = total_num_lesions++;
		total_num_lesions++;
    if (total_num_lesions>65000)
			err("total num lesions>65000");


    p[(border_size/2)*border_size+border_size/2]->set(border_size/2) ;
    //cells.push_back(c);
//    volume++ ;
    n=n0=1;
};

Lesion::~Lesion()
{
	for (unsigned int i=0; i < border_size*border_size; i++){
		delete p[i] ;
	};
	delete [] p ;
};


void Lesion::update_border_size()
{
  //int i,j,k;
  unsigned int new_border_size = int(border_size*1.25) ;
  if (new_border_size%2==1) new_border_size++ ; // make sure it's even
  unsigned int diff = int((new_border_size-border_size)/2) ;

#ifdef PUSHING
// fill in -1 in new sites  and extant lattice
  for (unsigned int i=0;i<border_size;i++) {
    for (unsigned int j=0;j<border_size;j++) {
      Sites *np=new Sites[new_border_size] ;
      for (unsigned int k=0;k<new_border_size;k++) np[k]=-1 ;
      for (unsigned int k=diff;k<border_size+diff;k++) np[k]=p[i][j][k-diff] ;
      delete p[i][j] ;
      p[i][j]=np ;
    }
  }

  Sites ***np=new Sites**[new_border_size] ;
  for (unsigned int i=0;i<new_border_size;i++) np[i]=new Sites*[new_border_size] ;

  for (unsigned int i=0;i<new_border_size;i++) {
    for (unsigned int j=0;j<new_border_size;j++) {
      if (i<diff || i>=border_size+diff || j<diff || j>=border_size+diff) {
        np[i][j]=new Sites[new_border_size] ;
        for (unsigned int k=0;k<new_border_size;k++) np[i][j][k]=-1 ;
      } else {
        np[i][j]=p[i-diff][j-diff] ;
      }
    }
  }

  for (unsigned int i=0;i<border_size;i++) delete [] p[i] ;

#else
// extant lattice
  for (unsigned int i=0;i<border_size*border_size;i++) {
    Sites *np=new Sites(new_border_size) ;
    for (unsigned int k=diff;k<border_size+diff;k++) if (p[i]->is_set(k-diff)) np->set(k) ;
    delete p[i] ;
    p[i]=np ;
  }
	if (float(new_border_size)*float(new_border_size)>2e9)
		err("new border size too large",new_border_size) ;
  Sites **np=new Sites*[new_border_size*new_border_size] ;

  for (unsigned int i=0;i<new_border_size;i++) {
    for (unsigned int j=0;j<new_border_size;j++) {
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

void Lesion::one_move_step(std::vector<Lesion*> lesions) {
  /*int i,j;*/
  double mthis = this->n ;
  for (unsigned int i=0;i<closest.size();i++) {
    vecd dr=lesions[closest[i]]->r - this->r ;
    double r2=squared(dr), sumrad2 = SQR(this->rad + lesions[closest[i]]->rad) ;
    if (r2 < sumrad2) {
      double mi=lesions[closest[i]]->n ;
      double disp = (sqrt(sumrad2 / r2) - 1) ;
//      if (fabs(disp) > maxdisp) maxdisp = fabs(disp);
//      if (fabs(disp) > maxdisp) maxdisp = fabs(disp);
      dr*=disp*1.1 ;
      this->r-=dr*mi/(mi+mthis) ;
      lesions[closest[i]]->r+=dr*mthis/(mi+mthis) ;
    }
  }
}

void Lesion::find_closest(std::vector<Lesion*> lesions)
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

//void Lesion::reduce_overlap(std::vector<Lesion*> lesions)
//{
//  int i,k,temp;
//  int *ind=new int[lesions.size()] ;
//  for (unsigned int j=0;j<lesions.size();j++) ind[j]=j ;
//  do {
//    maxdisp=0 ;
//    for (unsigned int j=0;j<lesions.size();j++)
//			{ k=_drand48()*lesions.size() ; SWAP(ind[j],ind[k]) ; }
//    for (unsigned int j=0;j<lesions.size();j++) {  // go through a random permutation
//      i=ind[j] ;
//      lesions[i]->one_move_step(lesions);
//
//      vecd dr=lesions[i]->r - lesions[i]->r_old ;
//      if (squared(dr)>SQR(lesions[i]->rad)) lesions[i]->find_closest(lesions);
//    }
//  } while (maxdisp>1e-2) ;
//  delete [] ind ;
//}

void cleanSimData(SimData & sim_data)
{
//std::cout<<"Hi1!!!"<<'\n';

//lesions preparation
	for (auto i:sim_data.lesions)
		if (i!=nullptr) delete i;
	sim_data.lesions.clear();
//genotypes preparation
	for (auto &i: sim_data.genotypes)
		delete i;
//if (i!=NULL)
//std::cout<<"Hi3!!!"<<'\n';
	sim_data.genotypes.clear();
//std::cout<<"Hi4!!!"<<'\n';
//cell preparation
	sim_data.cells.clear();

};

void initSimData(const float driver_adv, const float driver_mut_rate,
								 SimData & sim_data)
{
	cleanSimData(sim_data);

	//const float birth0 = 1;
	const float death0  = 1 - driver_adv;
  //max_growth_rate = growth0;
	//sim_data.genotypes.push_back(new Genotype(growth0, death0));

// cells preparation
  Cell cell;
  cell.x = cell.y = cell.z = 0;
//  cell.gen = sim_data.genotypes[0];
	//sim_data.root_genotype = new Genotype(growth0, death0);

	sim_data.genotypes.push_back(new Genotype(/*birth0,*/ death0));

// init lesions + first cell
  int total_num_lesions = 0;

  sim_data.lesions.push_back(new Lesion(cell, total_num_lesions));
  cell.gen = sim_data.genotypes[0];

  sim_data.cells.push_back(cell);
}

#ifdef MOORE_NEIGHBOURHOOD
const int _nonn=26 ;
const int kx[27]={0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1},
          ky[27]={0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1},
          kz[27]={0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1};
int kln[27] ; // this is filled with lengths of (kx,ky,kz)
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
  }
 // start_clock = clock();
}
//#endif
//Mykola end

//void end() {
//  fclose(times) ;
//}

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
  y=(border_size+y+ky[n])%border_size ;
  x=(border_size+x+kx[n])%border_size;
}


inline int free_sites(unsigned int n, std::vector<Cell> cells, Lesion *ll)
{
//  Lesion *ll=lesions[cells[n].lesion] ;
  int border_size = ll->border_size;
  int k = cells[n].x + border_size/2,
			j = cells[n].y + border_size/2,
			i = cells[n].z + border_size/2;
	return ll->no_free_sites(k,j,i);
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


//-----------------------------------------------------
//#if defined(NORMAL)

void checkEvoParameters(const float driver_adv,
										    const float mutation_rate,
										    const float driver_mutational_rate,
										    const unsigned int exit_size)
{
	if(	(0 > driver_adv) || (1 < driver_adv)	)
		err("a problem with the driver_adv parameter ", driver_adv);
	if(0 > mutation_rate)
		err("a problem with the mutation_rate parameter ", mutation_rate);
	if(	(0 > driver_mutational_rate) || (1	<	driver_mutational_rate)	)
		err("a problem with the driver_mutational_rate parameter ", driver_mutational_rate);
	if(0	>	exit_size)
		err("a problem with the exit_size parameter", exit_size);
};

inline void deleteGenotype(Genotype * gen,
													 std::vector<Genotype *> & genotypes)
{
	if (gen != genotypes.back()) {
		unsigned int index = gen->gen_array_index;
		genotypes[index] =  genotypes.back();
		genotypes[index] -> gen_array_index = index;
	};
	genotypes.pop_back();
	delete gen;
};

void cleanGenInfoFromVec(Genotype * current_gen,
												 std::vector<Genotype *> & genotypes)
{
	if (current_gen->getCellNum() == 0){
		if (current_gen->descendant_number == 0){
			current_gen->delMotherConnection();
			Genotype * ancestor_gen = current_gen->ancestor_gen;
			deleteGenotype(current_gen, genotypes);
			if (ancestor_gen != nullptr){
				cleanGenInfoFromVec(ancestor_gen, genotypes);
			};
		};
	};
};

int mainProgWellMixed(Evolution & evolution, SimData & sim_data)
{
	const float driver_adv						=	evolution.driver_adv;
	const float mutation_rate					=	evolution.mutation_rate;
	const float driver_mutation_rate	=	evolution.driver_mutation_rate;
	const unsigned int exit_size			=	evolution.stop_time_cell_num;
	checkEvoParameters(driver_adv, mutation_rate,
										 driver_mutation_rate, exit_size);
  unsigned int n; // index of a randomly chosen cell
  //double time = 0; // timer of the evolution

  // Total number of single nucleotide polymorphisms.
  // Can be increased during a sell devision (creation)
  // within new Genotype initiation: genotypes.push_back
	int total_num_snps = 0;

	const float max_death_rate = 1 - driver_adv;
	const float max_birth_rate = 1;
  for(;;) {
		int cell_vec_size = sim_data.cells.size();
		// Timer information. In cancer research timing of events is obsuted this simulations we don't know time
    n		 =	_drand48() * cell_vec_size;
    Genotype * current_gen = sim_data.cells[n].gen;


		float q = _drand48() * (max_death_rate + max_birth_rate);

		if (q < max_birth_rate) {
			Cell c;
			c.x	=	0;
			c.y	=	0;
			c.z	=	0;
			unsigned int num_snps = poisson(mutation_rate); // mutations for newly produced cell
			if (num_snps > 0) {// if number of SNP is not 0, then create a new genotype
				Genotype * new_genotype =	new Genotype(current_gen,
																	 sim_data.genotypes.size(),
																	 driver_adv,
																	 driver_mutation_rate,
																	 num_snps,
																	 total_num_snps);
				c.gen = new_genotype;
				sim_data.genotypes.push_back(new_genotype); // add new genotype
			} else { // if number of SNP is 0 add 1 cell in a correspondent genotype
				c.gen = current_gen;
				current_gen->incCellNum();
			};

			sim_data.cells.push_back(c); // add the new cell

			num_snps = poisson(mutation_rate); // SNPs old cell mutations
			if (num_snps > 0) {
				current_gen->decCellNum();
				Genotype * new_genotype =	new Genotype(current_gen,
																						 sim_data.genotypes.size(),
																						 driver_adv,
																						 driver_mutation_rate,
																						 num_snps,
																						 total_num_snps);
				sim_data.cells[n].gen = new_genotype;
				sim_data.genotypes.push_back(new_genotype); // add new genotype
				cleanGenInfoFromVec(current_gen, sim_data.genotypes);
			};
		} else {
//				if (_drand48() < current_gen->getDeath()){
			if((q>=1) && (q<(current_gen->getDeath() + max_birth_rate))){
				current_gen->decCellNum();
				cleanGenInfoFromVec(current_gen, sim_data.genotypes);

				if (n != (sim_data.cells.size()-1)){
					sim_data.cells[n]=sim_data.cells[sim_data.cells.size()-1];
				};
				sim_data.cells.pop_back();
			};
//				} else {
		};

		if (sim_data.cells.size() == 0) return 1;
		if (sim_data.cells.size() >= exit_size) {
				std::cout<<exit_size<<"; "<<sim_data.genotypes.size()<<"; "<</*genotype_death<<*/
				"; "<<driver_adv<<"; "<<mutation_rate <<"; !!!"<< driver_mutation_rate<<'\n';
				return 4 ;
		};
	};
}



int mainProg3D(Evolution & evolution, SimData & sim_data)
{
	const float driver_adv 						= evolution.driver_adv;
	const float mutation_rate 				= evolution.mutation_rate;
	const float driver_mutation_rate	= evolution.driver_mutation_rate;
	const unsigned int exit_size 			= evolution.stop_time_cell_num;

	checkEvoParameters(driver_adv,
										 mutation_rate,
										 driver_mutation_rate,
										 exit_size);

  unsigned int n; // index of a randomly chosen cell
  //double time = 0; // timer of the evolution
	//	short int total_number_lesions = 1;
	//	unsigned int genotype_death = 1;
  // Total number of single nucleotide polymorphisms.
  // Can be increased during a sell devision (creation)
  // within new Genotype initiation: genotypes.push_back
  int total_num_snps = 0;
	const float max_death_rate = 1 - driver_adv;
	const float max_birth_rate = 1;

////////////////////////////////////////////////////////////////////////////////
////description
////////////////////////////////////////////////////////////////////////////////
////  main_loop{
////  	time is change (according to tsc times number of cells ) and balance in step of birth event
////    take a cell randomly (n)
////    birth event is taken place and
////    {
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
	Lesion * ll = sim_data.lesions[0];
  for(;;) {      // main loop
		int cell_vec_size = sim_data.cells.size();
		if (cell_vec_size % 100000 == 0)
			std::cout << "cell_vec_size=" <<cell_vec_size<<'\n';
//    We don't need time
//    time += timescale / cell_vec_size;
    n=_drand48() * cell_vec_size; // choose number for a random current cell
    Genotype * current_gen = sim_data.cells[n].gen;
    //size of grid (wx - from -wx/2 to wx/2)
    unsigned int border_size = ll->border_size;
    int k = sim_data.cells[n].x + (border_size / 2);
    int j = sim_data.cells[n].y + (border_size / 2);
    int i = sim_data.cells[n].z + (border_size / 2);
    int need_border_size_update = 0;
		//if grid is too small change it
    if (k < 2 || k >= ((int)border_size - 3) ||
				j < 2 || j >= ((int)border_size - 3) ||
				i < 2 || i >= ((int)border_size - 3)) need_border_size_update = 1;
    if (ll->p[i * border_size + j]->is_set(k) == 0)
			err("ll->p[i][j][k]==0, border size =", border_size);

		float q = _drand48() * (max_death_rate + max_birth_rate);
		float br;

    if (free_sites(n, sim_data.cells, ll)>0)
			br = 1;
		else
			br = 0;

		float dr = sim_data.cells[n].gen->getDeath();

//create new cell
		if (q < br) {

      int kn=k, jn=j, in=i ;
      ll->choose_nn(kn,jn,in) ;

			Cell c;
			c.x = kn - (border_size / 2);
			c.y = jn - (border_size / 2);
			c.z = in - (border_size / 2);
			//c.lesion = sim_data.cells[n].lesion;
			ll -> p[in * border_size + jn]->set(kn) ; // set grid occupied
			int num_snps = poisson(mutation_rate) ; // newly produced cell mutants
			if (num_snps > 0) {// num of SNP is not 0 => create and add new genotype
				Genotype * new_genotype =	new Genotype(current_gen,
																	 sim_data.genotypes.size(),
																	 driver_adv,
																	 driver_mutation_rate,
																	 num_snps,
																	 total_num_snps);
				c.gen = new_genotype;
				sim_data.genotypes.push_back(new_genotype); // add new genotype
			} else { // if number of SNP is 0 add 1 cell in a correspondent genotype
				c.gen = current_gen;
				current_gen->incCellNum();
			};
			sim_data.cells.push_back(c); // add the newly created cell
			ll->n++ ;

			num_snps = poisson(mutation_rate); // SNPs old cell mutations
			if (num_snps > 0) {
				current_gen->decCellNum();
				Genotype * new_genotype =	new Genotype(current_gen,
																						 sim_data.genotypes.size(),
																						 driver_adv,
																						 driver_mutation_rate,
																						 num_snps,
																						 total_num_snps);
				sim_data.cells[n].gen = new_genotype;
				sim_data.genotypes.push_back(new_genotype); // add new genotype
				cleanGenInfoFromVec(current_gen, sim_data.genotypes);
			};
		} else {
			if((q>=br) && (q<(br + dr))){
// now we implement death event
				if(_drand48() < current_gen->getDeath()){
					ll->p[ i * border_size+j] -> unset(k);
					ll->n--;
					current_gen->decCellNum();
					//std::cout<<"1to"<<'\n';
					cleanGenInfoFromVec(current_gen, sim_data.genotypes);
					if (n != (sim_data.cells.size()-1)) {
						sim_data.cells[n]=sim_data.cells[sim_data.cells.size()-1];
					};
					sim_data.cells.pop_back();
				} else {
					// chosen cell is survived and a mutation occur
					int num_snps = poisson(mutation_rate) ; // old cell mutates
					if (num_snps > 0) {
						current_gen->decCellNum();
						Genotype * new_genotype =	new Genotype(current_gen,
																									 sim_data.genotypes.size(),
																									 driver_adv,
																									 driver_mutation_rate,
																									 num_snps,
																									 total_num_snps);
						sim_data.cells[n].gen = new_genotype;
						sim_data.genotypes.push_back(new_genotype); // add new genotype
						cleanGenInfoFromVec(current_gen, sim_data.genotypes);
					};
					//chosenCellSurviveAndMutate(current_gen, sim_data);
				};
			};
		};

		if (need_border_size_update && ll!=NULL) ll->update_border_size() ;

		if (sim_data.cells.size()==0) return 1 ;
		if (sim_data.cells.size()>=exit_size) {
			std::cout<<exit_size<<"; "<<sim_data.genotypes.size()<<"; "<<
			"; "<<driver_adv<<"; "<<mutation_rate <<"; "<< driver_mutation_rate<<'\n';
			return 4 ;
		};
	};
};
//#endif // NORMAL


//-----------------------------------------------------------------------------

#if defined(FASTER_KMC) || defined(GILLESPIE)

int mainProg3D(int exit_size, int save_size, double max_time, double wait_time)
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
    double tot_rate=cells.size()*(max_growth_rate + max_death_rate) ;
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
#ifdef PUSHNG
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
    if (ntot>=exit_size) return 4 ;

  }

}

#endif // FASTER_KMC or GILLESPIE
