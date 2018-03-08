#ifndef PROG_PARAMETERS_H_INCLUDED
#define PROG_PARAMETERS_H_INCLUDED

//external parameters (from configuration files)
#include <conf_structures.h>
#include <string>
#include <iostream>
//#include <random>


struct ParPoint{
  float driver_adv;
  float mut_rate;  // gama; // these are rates per daughter cell.
                   //Rates per diploid exome will be 2x higher
                   //(these values are given in the paper)
  float driver_mut_rate; //driver_prob;
};


struct Results {
	float tree_dist;
	float dist;
	float balk_dist;
	float weight;

//	int delta_mut;
//	int num_mut;
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
void SimulateAndWrite2File(const ParPoint par_point,
													 const AllSimParameters & conf_settings,
													 const int rand);

int GetNumLinesInComparisonFile(const AllSimParameters & conf_settings);


//void GetTreeComparisonVecFromFile(const Parameters & pars,
//													        std::vector<TreeComparison> & tree_comparison_vec);
void GetTreeComparisonVecFromFile(const AllSimParameters & conf_settings, const int diff,
																	std::vector<TreeComparison> & tree_comparison_vec);
void SetTreeComparisonVec2File(const std::vector<TreeComparison> & tree_comparison_vec,
															 const AllSimParameters & conf_settings, const int number_of_iteration);
void AddToLogFile(const AllSimParameters & conf_settings, const unsigned int iteration_num,
									const unsigned int trials, const TreeComparison & relative_error,
									const TreeComparison & variance);

std::string GetPathToEtalonResultCatalogue(const AllSimParameters & conf_settings);

std::string GetFullComparisonFileName(const AllSimParameters & conf_settings);
std::string GetFullProgName();
int GetSearchParNumber();
int GetRandNumber(Vars & vars);
float GetLeftBorderDriverAdvantage(const AllSimParameters & conf_settings);
float GetRigthBorderDriverAdvantage(const AllSimParameters & conf_settings);
float GetLeftBorderMutationRate(const AllSimParameters & conf_settings);
float GetRigthBorderMutationRate(const AllSimParameters & conf_settings);
float GetLeftBorderDriverMutationRate(const AllSimParameters & conf_settings);
float GetRigthBorderDriverMutationRate(const AllSimParameters & conf_settings);

#endif // PROG_PARAMETERS_H_INCLUDED
