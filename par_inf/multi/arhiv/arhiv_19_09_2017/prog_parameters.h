#ifndef PROG_PARAMETERS_H_INCLUDED
#define PROG_PARAMETERS_H_INCLUDED

//external parameters (from configuration files)
#include <conf_structures.h>
//#include <random>
//#include <string>


struct ParPoint{
  float driver_adv;
  float mut_rate;  // gama; // these are rates per daughter cell.
                   //Rates per diploid exome will be 2x higher
                   //(these values are given in the paper)
  float driver_mut_rate; //driver_prob;
};


struct Results {
	float tree_dist;
	int delta_mut;
	int num_mut;
	float weight;
};

struct TreeComparison{
	ParPoint par_point;
	Results results;
	TreeComparison & operator = (const TreeComparison& a) {
		par_point = a.par_point;
    results = a.results;
		return *this;
	}
};

struct Algo {
    int parameter_population_num;
    int parameter_population_num_alpha;
    float p_acc_min;
    float final_epsilon;
    float final_iteration;
};

struct Vars {
	std::mt19937 gen;
	float max_dist;
	float p_acc;
	unsigned n_trials;
};

struct InfoSearch {
	float min_;
	float max_;
  float true_val;
  bool is_mutable;
};

struct SearchSpace{
  InfoSearch driver_adv;
  InfoSearch mut_rate;
  InfoSearch driver_mut_rate;
};

struct Parameters {
//	string catalog_file_name;
//	string par_file_name;
//	string comparison_file_name;
//	string log_file_name;

//	int seed; //for initialization of Mersenne Twister pseudo-random number gen
//  SearchSpace search_space;
//	int number_of_parameters;
//  Algo algo;
	AllSimParameters conf_settings;
  Vars vars;
};

//augment functions

bool CompareTwoTrees(const TreeComparison tree_results_1,
										 const TreeComparison tree_results_2);


// the file is named by "%s/file_tree_comparison.dat" from cancer_3d
void SimulateAndWrite2File(const ParPoint par_point, Parameters & pars);

int GetNumLinesInComparisonFile(Parameters & pars);


//void GetTreeComparisonVecFromFile(const Parameters & pars,
//													        std::vector<TreeComparison> & tree_comparison_vec);
void GetTreeComparisonVecFromFile(const Parameters & pars, const int diff,
																	std::vector<TreeComparison> & tree_comparison_vec);
void SetTreeComparisonVec2File(const std::vector<TreeComparison> & tree_comparison_vec,
															 const Parameters & pars, const int number_of_iteration);
void AddToLogFile(const Parameters & pars,const unsigned int iteration_num,
									const unsigned int trials, const TreeComparison & relative_error,
									const TreeComparison & variance);

std::string GetPathToEtalonResultCatalogue(const Parameters & pars);

std::string GetFullComparisonFileName(Parameters & pars);

int GetSearchParNumber();

float GetLeftBorderDriverAdvantage(const Parameters & pars);
float GetRigthBorderDriverAdvantage(const Parameters & pars);
float GetLeftBorderMutationRate(const Parameters & pars);
float GetRigthBorderMutationRate(const Parameters & pars);
float GetLeftBorderDriverMutationRate(const Parameters & pars);
float GetRigthBorderDriverMutationRate(const Parameters & pars);

#endif // PROG_PARAMETERS_H_INCLUDED
