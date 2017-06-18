
//i/o includes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//#include <boost/filesystem.hpp>

//parallel calculations
#include <omp.h>

#include <random>
#include <algorithm>    // std::sort
#include <vector>

//#include <ctime>
//#include "stdafx.h"
#include <fstream>

#include <boost/limits.hpp>
// /multivariate_gaussian_distribution/ //
#include "multivariate_gaussian_distribution.h"
#include "cholesky.hpp"
#include "inverse_matrix.h"

#include "debug_functions.h"
//#define DEBUG_REGIME  // check regime





namespace ublas = boost::numeric::ublas;

using namespace std;
//---------------conf file--------------------------------
//begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
//end
//---------------conf file--------------------------------


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

struct TreeComparison {
	ParPoint par_point;
	Results results;
	TreeComparison & operator = (const TreeComparison& a) {
		par_point = a.par_point;
    results =a.results;
		return *this;
	}
};

struct Algo {
    int parameter_population_num;
    int parameter_population_num_alpha;
    float p_acc_min;
    float final_epsilon;
};

struct Vars {
	mt19937 gen;
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
	string catalog_file_name;
	string par_file_name;
	string comparison_file_name;
	int seed; //for initialization of Mersenne Twister pseudo-random number gen
  SearchSpace search_space;
	int number_of_parameters;
  Algo algo;
  Vars vars;
};

//augment functions
void err(const char *reason); // exits program with correspondent reason

bool comparison_of_two_trees(TreeComparison tree_results_1,
														 TreeComparison tree_results_2) {
    return (tree_results_1.results.tree_dist <
						  tree_results_2.results.tree_dist);
};

// the file is named by "%s/file_tree_comparison.dat" from cancer_3d
void SimulateAndWrite2File(const ParPoint par_point, Parameters & pars) {
    char command_line[256];
		std::random_device rd;
		std::uniform_int_distribution<int> rnd(1,2000000000);

    int rand = rnd(pars.vars.gen);
#if defined __linux
    sprintf(command_line,"./cancer_3d.exe %s %f %f %f %d",
						pars.par_file_name.c_str(),
						par_point.driver_adv,
						par_point.mut_rate,
						par_point.driver_mut_rate,
						rand);
#else
    sprintf(command_line,"cancer_3d.exe %s %f %f %f %d",
						pars.par_file_name.c_str(),
						par_point.driver_adv,
						par_point.mut_rate,
						par_point.driver_mut_rate,
						rand);
#endif
    //cout <<"command_line[256]="<<command_line<<endl;
    system (command_line);
};

int get_num_lines_in_comp_file(Parameters & pars)
{
    std::ifstream file_tree_comparison;
    char full_name_file[256];

    sprintf(full_name_file,"%s/%s",pars.catalog_file_name.c_str(),
															pars.comparison_file_name.c_str());
    file_tree_comparison.open(full_name_file);

    int number_of_lines = 0;
    std::string line;
    while (std::getline(file_tree_comparison, line)) number_of_lines++;
    return number_of_lines;
};

void GetTreeComparisonVecFromFile(std::vector<TreeComparison> &
																	tree_comparison_vec,
														      Parameters & pars) {
	std::ifstream file_tree_comparison;
	char name_file[256];
	sprintf(name_file,"%s/%s", pars.catalog_file_name.c_str(),
					                   pars.comparison_file_name.c_str());
	file_tree_comparison.open(name_file);
	if (file_tree_comparison.is_open()==false) {
		cout<<"err..."<<endl;
	} else {
		for(int i = 0; i < pars.algo.parameter_population_num; i++) {
			TreeComparison A;
			file_tree_comparison >> A.par_point.driver_adv
													 >> A.par_point.mut_rate
													 >> A.par_point.driver_mut_rate
													 >> A.results.tree_dist
													 >> A.results.delta_mut
													 >> A.results.num_mut;
			if (!file_tree_comparison.eof()) tree_comparison_vec.push_back(A);
		};
 };
}

void GetTreeComparisonVecFromFile(std::vector<TreeComparison> &
																	tree_comparison_vec,
														      Parameters & pars, int diff) {
	std::ifstream file_tree_comparison;
	char name_file[256];
	sprintf(name_file,"%s/%s", pars.catalog_file_name.c_str(),
					                   pars.comparison_file_name.c_str());
	file_tree_comparison.open(name_file);
	if (file_tree_comparison.is_open()==false) {
		cout<<"err..."<<endl;
	} else {
		for(int i = 0; i < diff; i++) {
			TreeComparison A;
			file_tree_comparison >> A.par_point.driver_adv
													 >> A.par_point.mut_rate
													 >> A.par_point.driver_mut_rate
													 >> A.results.tree_dist
													 >> A.results.delta_mut
													 >> A.results.num_mut;
			if (!file_tree_comparison.eof()) tree_comparison_vec.push_back(A);
		};
 };
}

void InitWeight(std::vector <TreeComparison> & tree_comparison_vec,
					 int parameter_population_num) {
	for(int i=0; i < parameter_population_num; i++) {
		tree_comparison_vec[i].results.weight = 1.0/parameter_population_num;
	};
};

void SetTreeComparisonVec2File(std::vector<TreeComparison> &
															 tree_comparison_vec,
															 Parameters & pars,
															 int number_of_iteration) {
	std::ofstream file_tree_comparison;
	char name_file[256];

	sprintf(name_file,"%s/%d_%s", pars.catalog_file_name.c_str(),
																number_of_iteration,
																pars.comparison_file_name.c_str());
	file_tree_comparison.open(name_file);

	if (file_tree_comparison.is_open()==false) cout<<"err..."<<endl;

	file_tree_comparison.precision(5);
	file_tree_comparison.setf(std::ios::fixed, std:: ios::floatfield);

	for(unsigned int i=0; i < tree_comparison_vec.size(); i++)
		file_tree_comparison<<
			tree_comparison_vec[i].par_point.driver_adv<<"\t"<<
			tree_comparison_vec[i].par_point.mut_rate<<"\t"<<
			tree_comparison_vec[i].par_point.driver_mut_rate<<"\t"<<
			tree_comparison_vec[i].results.tree_dist<<"\t"<<
			tree_comparison_vec[i].results.delta_mut<<"\t"<<
			tree_comparison_vec[i].results.num_mut<<"\t"<<
			tree_comparison_vec[i].results.weight<<endl;
};

void InitLoop(Parameters & pars) {
	char name_file[256];
	/*file_tree_comparison.dat*/
	sprintf(name_file,"%s/%s", pars.catalog_file_name.c_str(),
					                   pars.comparison_file_name.c_str());
	remove(name_file);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> driver_adv_gen(
																				pars.search_space.driver_adv.min_,
																				pars.search_space.driver_adv.max_ );
	std::uniform_real_distribution<float> mut_rate_gen(
																				pars.search_space.mut_rate.min_,
																				pars.search_space.mut_rate.max_);
	std::uniform_real_distribution<float> driver_mut_rate_gen(
																				pars.search_space.driver_mut_rate.min_,
																				pars.search_space.driver_mut_rate.max_);

	do {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  		  #pragma omp parallel for
			  for(int i = 0; i <= pars.algo.parameter_population_num; i++) {
					//ParPoint random_par_point = ;
					SimulateAndWrite2File({driver_adv_gen(pars.vars.gen),
																 mut_rate_gen(pars.vars.gen),
																 driver_mut_rate_gen(pars.vars.gen)}, pars);
			  };
	} while (get_num_lines_in_comp_file(pars) <=
					 pars.algo.parameter_population_num);

};

//parameter_population_num_alpha
void WeightNormalization (std::vector<TreeComparison>& tree_comparison_vec,
													 int begin_index, int end_index) {
	float weight_sum_alpha = 0;
	for(int i = begin_index; i<end_index; i++)
		weight_sum_alpha += tree_comparison_vec[i].results.weight;

 // cout<< "parameter_population_num_alpha="<<parameter_population_num_alpha;
 // cout<< "weight_sum_alpha="<<weight_sum_alpha;
	for(int i=begin_index; i<end_index; i++)
		tree_comparison_vec[i].results.weight=
			((tree_comparison_vec[i].results.weight*1.0) / (1.0*weight_sum_alpha));
};

ublas::vector<float> op(const ublas::vector<float> & weighted_vector,
												TreeComparison T) {
    ublas::vector<float> w_v(weighted_vector.size());
    w_v(0) = weighted_vector(0) + T.results.weight * T.par_point.driver_adv;
    w_v(1) = weighted_vector(1) + T.results.weight * T.par_point.mut_rate;
    w_v(2) = weighted_vector(2) + T.results.weight * T.par_point.driver_mut_rate;
    return w_v;
};

float SumOfSquareWeights(std::vector<TreeComparison> & array_of_tree_comparisons,
												 int parameter_population_num_alpha) {
  float sum = 0;
  for (int i = 0; i<parameter_population_num_alpha; i++)
		sum += array_of_tree_comparisons[i].results.weight*
					   array_of_tree_comparisons[i].results.weight;
  return sum;
};

ublas::vector<float> GetWeightedMeanVec(
			 std::vector<TreeComparison> & tree_comparison_vec, Parameters & pars) {
  int size_of_matrix = pars.number_of_parameters;
  int parameter_population_num_alpha = pars.algo.parameter_population_num_alpha;
  ublas::vector<float> weighted_mean_vector(size_of_matrix);

  weighted_mean_vector = ublas::zero_vector<float> (size_of_matrix);
  weighted_mean_vector = std::accumulate(
                            tree_comparison_vec.begin(),
                            tree_comparison_vec.begin() +
																				parameter_population_num_alpha,
                            weighted_mean_vector,
					                  op);
	if (!(pars.search_space.driver_adv.is_mutable)) weighted_mean_vector[0] =
		pars.search_space.driver_adv.true_val;
	if (!(pars.search_space.mut_rate.is_mutable)) weighted_mean_vector[1] =
		pars.search_space.mut_rate.true_val;
	if (!(pars.search_space.driver_mut_rate.is_mutable)) weighted_mean_vector[2] =
		pars.search_space.driver_mut_rate.true_val;

  return weighted_mean_vector;

};


ublas::matrix<float, ublas::row_major> GetTwiceWeightedEmpiricalCovariance(
												ublas::vector<float> & weighted_mean_vector,
												std::vector<TreeComparison> & array_of_tree_comparisons,
												Parameters & pars) {
  int size_matrix = pars.number_of_parameters;
  int parameter_population_num_alpha = pars.algo.parameter_population_num_alpha;
  //init epsion_matrix
  ublas::matrix<float, ublas::row_major> epsilon_matrix(size_matrix, size_matrix);
  epsilon_matrix = ublas::zero_matrix<float>(size_matrix, size_matrix);
#if defined(DEBUG_REGIME)
    std::cout<<"In function GetTwiceWeightedEmpiricalCovariance:"<<std::endl;
#endif

  for (int i=0; i<parameter_population_num_alpha; i++)	{
    ublas::vector<float> parameter_vector(size_matrix);

    parameter_vector(0) = array_of_tree_comparisons[i].par_point.driver_adv;
    parameter_vector(1) = array_of_tree_comparisons[i].par_point.mut_rate;
    parameter_vector(2) = array_of_tree_comparisons[i].par_point.driver_mut_rate;

    ublas::vector<float> diff_mean_vector(size_matrix);

    diff_mean_vector[0] = parameter_vector[0] - weighted_mean_vector[0];
    diff_mean_vector[1] = parameter_vector[1] - weighted_mean_vector[1];
    diff_mean_vector[2] = parameter_vector[2] - weighted_mean_vector[2];
    //for Choletski decomposition (create positive definite matrices)
    //if (!(pars.search_space.driver_adv.is_mutable)) diff_mean_vector[0] = 3e-6;
    //if (!(pars.search_space.mut_rate.is_mutable)) diff_mean_vector[1] = 2e-6;
    //if (!(pars.search_space.driver_mut_rate.is_mutable)) diff_mean_vector[2] = 1e-6;
    epsilon_matrix += array_of_tree_comparisons[i].results.weight *
                       (ublas::outer_prod(diff_mean_vector,diff_mean_vector));

  };

  float sum_sq_wei = SumOfSquareWeights(array_of_tree_comparisons,
																					 parameter_population_num_alpha);
  if (sum_sq_wei<1) {
		epsilon_matrix*= ( 2.0 / (1 - sum_sq_wei));
	} else {
		err("in func GetTwiceWeightedEmpiricalCovariance sum_sq_wei>=1");
	};

 //for Choletski decomposition (create positive definite matrices)
  for (unsigned int i=0;i<epsilon_matrix.size1();i++)
		if (epsilon_matrix(i,i)==0) epsilon_matrix(i,i) = 1e-12;

  for(unsigned int i = 0; i < epsilon_matrix.size1(); i++)
		for(unsigned int j = 0; j < epsilon_matrix.size2(); j++)
				if (epsilon_matrix(i,j)==0) epsilon_matrix(i,j) = 1e-15;

  return epsilon_matrix;
};

int get_index_of_parameter_vector (
		std::vector<TreeComparison> array_of_tree_comparisons,
		int parameter_population_num_alpha, mt19937 & gen) {
	std::uniform_real_distribution<float> dis(0, 1);

	float uniform_coin=dis(gen);
	float sum = array_of_tree_comparisons[0].results.weight;
	int indicator = 0;
	while (sum < uniform_coin){
		sum += array_of_tree_comparisons[indicator].results.weight;
		indicator++;
	};
	return indicator;
};


float get_reflected_coord (float old_coord, float l_b/*left border*/,
																						float r_b/*left border*/) {
	//std::cout<<"old_coord="<<old_coord<<"left_border="<<l_b<<"right_border="<<r_b<<endl;
	float diff_min = old_coord - l_b;
	if (diff_min<0) {
		return (l_b + abs(diff_min));
	} else {
		float diff_max = old_coord - r_b;
		if (diff_max<0) err("get_reflected_coord: diff_max<0 ");
		return  (r_b - diff_max);
	};

};

// reflect coords of parameters
// until all of them will be in acceptable parameter range
ublas::vector<float> ReflectedCoords (ublas::vector<float> & new_parameter_vector,
																			Parameters & pars) {
  std::vector <bool> is_in_dom_vec = {false, false, false};
  do {
		if ((pars.search_space.driver_adv.min_<=new_parameter_vector[0] &&
			pars.search_space.driver_adv.max_>=new_parameter_vector[0]) || is_in_dom_vec[0] ) {
			  is_in_dom_vec[0] = true;
		} else {
			//std::cout<<"ReflectedCoords driver_adv"<<std::endl;
			new_parameter_vector[0] = get_reflected_coord(new_parameter_vector[0],
															   pars.search_space.driver_adv.min_,
															   pars.search_space.driver_adv.max_);
		}; // process driver_adv coord
		if ((pars.search_space.mut_rate.min_<=new_parameter_vector[1] &&
				pars.search_space.mut_rate.max_>=new_parameter_vector[1]) || is_in_dom_vec[1]) {
			is_in_dom_vec[1] = true;
		} else {
			//std::cout<<"ReflectedCoords mut_rate"<<std::endl;
			new_parameter_vector[1] = get_reflected_coord(new_parameter_vector[1],
															   pars.search_space.mut_rate.min_,
															   pars.search_space.mut_rate.max_);
		}; // process mut_rate coord
		if ((pars.search_space.driver_mut_rate.min_<=new_parameter_vector[2] &&
				pars.search_space.driver_mut_rate.max_>=new_parameter_vector[2]) || is_in_dom_vec[2]) {
			is_in_dom_vec[2] = true;
		} else {
			//std::cout<<"ReflectedCoords driver_mut_rate"<<std::endl;
			new_parameter_vector[2] = get_reflected_coord(new_parameter_vector[2],
															   pars.search_space.driver_mut_rate.min_,
															   pars.search_space.driver_mut_rate.max_);
		}; // process driver_mut_rate coord
	} while (is_in_dom_vec[0] == false ||
					 is_in_dom_vec[1] == false ||
					 is_in_dom_vec[2] == false );
	return new_parameter_vector;
};
//    return parameter_vector;
ublas::vector<float> GetMultiNormalDistSample(
												ublas::matrix<float, ublas::row_major> & L,
                        ublas::vector<float> & parameter_vector,
                        Parameters & pars) {
	int size_matrix = L.size1();
	ublas::vector<float> new_parameter_vector(size_matrix);
	std::normal_distribution<float> n_d(0,1);

	new_parameter_vector[0] = n_d(pars.vars.gen);
	new_parameter_vector[1] = n_d(pars.vars.gen);
	new_parameter_vector[2] = n_d(pars.vars.gen);

	new_parameter_vector = prod(new_parameter_vector, L);
#if defined(DEBUG_REGIME)
	std::cout<<"We are in GetMultiNormalDistSample"<<std::endl;
  std::cout<<"before adding present param new_parameter_vector="<<
		new_parameter_vector<<std::endl;
#endif

	if (!(pars.search_space.driver_adv.is_mutable)) new_parameter_vector[0] = 0;
	if (!(pars.search_space.mut_rate.is_mutable)) new_parameter_vector[1] = 0;
  if (!(pars.search_space.driver_mut_rate.is_mutable)) new_parameter_vector[2] = 0;
	new_parameter_vector += parameter_vector;

	// reflect coords of parameters
	// until all of them will be in acceptable parameter range
	ReflectedCoords(new_parameter_vector, pars);
#if defined(DEBUG_REGIME)
	std::cout<<"After reflection new_parameter_vector="<<
		new_parameter_vector<<std::endl;
#endif

	return new_parameter_vector;
};

void SetProjs(ublas::vector<float> & proj_i_vector,
							ublas::vector<float> & proj_j_vector,
							std::vector<TreeComparison> & array_of_tree_comparisons,
							std::vector <bool> is_par_vec, int i, int j) {
	int index = 0;
	if (is_par_vec[0]) {
		proj_i_vector(index) = array_of_tree_comparisons[i].par_point.driver_adv;
		proj_j_vector(index) = array_of_tree_comparisons[j].par_point.driver_adv;
		index++;
	};
	if (is_par_vec[1]) {
		proj_i_vector(index) = array_of_tree_comparisons[i].par_point.mut_rate;
		proj_j_vector(index) = array_of_tree_comparisons[j].par_point.mut_rate;
		index++;
	};
	if (is_par_vec[2]) {
	 proj_i_vector(index) = array_of_tree_comparisons[i].par_point.driver_mut_rate;
	 proj_j_vector(index) = array_of_tree_comparisons[j].par_point.driver_mut_rate;
	};
};

float GetGaussianKernel(
				std::vector<TreeComparison> & array_of_tree_comparisons,
				ublas::matrix<float, ublas::row_major> & epsilon_matrix,
				int i, int j, Parameters & pars) {
#if defined(DEBUG_REGIME)
  std::cout<<"In GetGaussianKernel="<<std::endl;
#endif

	std::vector <bool> is_par_vec;
	is_par_vec.push_back(pars.search_space.driver_adv.is_mutable);
	is_par_vec.push_back(pars.search_space.mut_rate.is_mutable);
	is_par_vec.push_back(pars.search_space.driver_mut_rate.is_mutable);
	int num_par = 0;
	for (unsigned int it = 0; it < is_par_vec.size(); it++) if (is_par_vec[it]) num_par++;
	range r(0, num_par);
	ublas::matrix<float, ublas::row_major> proj_epsilon_matrix(num_par, num_par),
																	inverse_proj_epsilon_matrix(num_par, num_par);
	ublas::vector<float> proj_i_vector(num_par), proj_j_vector(num_par),
					             diff_proj_i_j_vector(num_par), result_of_mult(num_par);
  SetProjs(proj_i_vector,proj_j_vector,array_of_tree_comparisons,is_par_vec,i,j);

	diff_proj_i_j_vector = proj_i_vector - proj_j_vector;
#if defined(DEBUG_REGIME)
  std::cout<<"proj_i_vector="<<proj_i_vector<<std::endl;
  std::cout<<"proj_j_vector="<<proj_j_vector<<std::endl;
  std::cout<<"diff_proj_i_j_vector="<<diff_proj_i_j_vector<<std::endl;
#endif

  int index_row = 0;
  for(int row = 0; row < pars.number_of_parameters; row++) {
		if (is_par_vec[row]) {
			int index_col = 0;
			for(int col = 0; col < pars.number_of_parameters; col++) {
				if(is_par_vec[col]) {
					proj_epsilon_matrix(index_row, index_col) = epsilon_matrix(row, col);
					index_col++;
				};
			};
			index_row++;
		};
  };

	InvertMatrix(proj_epsilon_matrix, inverse_proj_epsilon_matrix);
	result_of_mult = prod(inverse_proj_epsilon_matrix, diff_proj_i_j_vector);
#if defined(DEBUG_REGIME)
	std::cout<<"epsilon_matrix=" << epsilon_matrix <<std::endl;
	std::cout<<"proj_epsilon_matrix" << proj_epsilon_matrix <<std::endl;
	std::cout<<"inverse_proj_epsilon_matrix="<<
						  inverse_proj_epsilon_matrix<<std::endl;
	std::cout<<"result_of_mult" << result_of_mult <<std::endl;
	std::cout<<"return"<< inner_prod(diff_proj_i_j_vector, result_of_mult)<<std::endl;
	int i_00;
	cin>>i_00;
#endif
	return inner_prod(diff_proj_i_j_vector, result_of_mult);
};

void SetNewWeights2NewTreeComparisons(
												std::vector<TreeComparison>& tree_comparison_vec,
												ublas::matrix<float, ublas::row_major> & epsilon_matrix,
												Parameters & pars,
												int index_of_begin) {
	int index_of_end = tree_comparison_vec.size();
	for (int i = index_of_begin; i < index_of_end; i++) {
		float sum_denomination = 0;
		for(int j=0; j<index_of_begin; j++) {
			sum_denomination += tree_comparison_vec[j].results.weight *
			GetGaussianKernel(tree_comparison_vec, epsilon_matrix, i, j, pars);
		};
		tree_comparison_vec[i].results.weight =
			1.0/(index_of_end * sum_denomination * 1.0);
//	  std::cout<<"we are in SetNewWeights2NewTreeComparisons" << std::endl <<
//	  "current index=" << i << "; index_of_end="<<index_of_end<<
//		"; sum_denomination=" << sum_denomination <<
//		";tree_comparison_vec[i].results.weight=" <<
//		tree_comparison_vec[i].results.weight << std::endl;
	};
};

void SubInnerLoop(Parameters & pars,
									std::vector<TreeComparison> & tree_comparison_vec,
									ublas::matrix<float, ublas::row_major> & epsilon_matrix,
									ublas::matrix<float, ublas::row_major> Cholesky_matrix) {
	int matrix_size = epsilon_matrix.size1();
	int index = get_index_of_parameter_vector(tree_comparison_vec,
																			pars.algo.parameter_population_num_alpha,
																			pars.vars.gen);
	ublas::vector<float> parameter_vector(matrix_size);
	parameter_vector[0] = tree_comparison_vec[index].par_point.driver_adv;
	parameter_vector[1] = tree_comparison_vec[index].par_point.mut_rate;
	parameter_vector[2] = tree_comparison_vec[index].par_point.driver_mut_rate;

	parameter_vector = GetMultiNormalDistSample(Cholesky_matrix,
																							parameter_vector,
																							pars);
  ParPoint par_point={ parameter_vector [0],
											 parameter_vector [1],
											 parameter_vector [2] };
	SimulateAndWrite2File(par_point, pars);
	pars.vars.n_trials++;
};

void InnerLoop(Parameters & pars,
							 std::vector<TreeComparison> & tree_comparison_vec,
							 ublas::matrix<float, ublas::row_major> & epsilon_matrix) {
	int matrix_size = epsilon_matrix.size1();
	pars.vars.max_dist = tree_comparison_vec[
												 tree_comparison_vec.size() -
	                       pars.algo.parameter_population_num_alpha - 1
																						].results.tree_dist;
  // init file for new loop
	char name_file[256];
	sprintf(name_file,"%s/%s",pars.catalog_file_name.c_str(),
					                  pars.comparison_file_name.c_str());
	remove(name_file);  // end init file for new loop
  //Cholesky_decompose
	ublas::matrix<float, ublas::row_major> L (matrix_size, matrix_size);
	L=ublas::zero_matrix<float>(matrix_size, matrix_size);
	cholesky_decompose(epsilon_matrix, L);
	// end Cholesky_decompose
#if defined(DEBUG_REGIME)
	cout<<"In InnerLoop"<<std::endl;
	cout<<"cholesky_decompose="<<L<<std::endl;
#endif

  int diff = pars.algo.parameter_population_num -
						 pars.algo.parameter_population_num_alpha;
	do{//shared(flag) #pragma omp parallel for
			int boundary = 2.0 * pars.algo.parameter_population_num;
  		//download new mingw in a case of problems for OpenMP
      #pragma omp parallel for
			for(int i = 0; i < boundary; i++) {
				SubInnerLoop(pars, tree_comparison_vec, epsilon_matrix, L);
			};
	} while(get_num_lines_in_comp_file(pars) <= diff);

  std::vector<TreeComparison> tree_comparison_vec_plus;
  GetTreeComparisonVecFromFile(tree_comparison_vec_plus, pars, diff);
	//send results to main array array_of_tree_comparisons from additional array
	for(int i = pars.algo.parameter_population_num_alpha;
			i < pars.algo.parameter_population_num;
			i++)
			tree_comparison_vec[i] = tree_comparison_vec_plus[i -
		                           pars.algo.parameter_population_num_alpha];

	SetNewWeights2NewTreeComparisons(tree_comparison_vec, epsilon_matrix, pars,
																		 pars.algo.parameter_population_num_alpha);
};




void MainLoop(Parameters & pars, std::vector<TreeComparison>& tree_comparison_vec);
Parameters InitMainParameters(int argc, char *argv[]);
void DoesInitiationPart(Parameters & pars);
void GetPars(Parameters & pars);


int main (int argc, char *argv[]) {
	Parameters pars = InitMainParameters(argc, argv);
	DoesInitiationPart(pars);
	std::vector<TreeComparison> tree_comparison_vec;
 GetTreeComparisonVecFromFile(tree_comparison_vec, pars);
 InitWeight(tree_comparison_vec, pars.algo.parameter_population_num);
#if defined(DEBUG_REGIME)
	SetTreeComparisonVec2File(tree_comparison_vec, pars, 100000);
#endif
  MainLoop(pars, tree_comparison_vec);
  return 0;
};

void err(const char *reason) {
  cout <<reason<<endl ;
#ifdef __WIN32
  system("pause") ;
#endif
  exit(0) ;
}

bool FileExists(const char *name) {
    ifstream my_file(name);
    if(my_file.fail()){ //File does not exist
        return false;
    } else { //otherwise, file exists
        return true;
    };
};


Parameters InitMainParameters(int argc, char *argv[]) {
  Parameters pars;
  const int expected_arg_num = 3; // the number of arguments in command line
  if (argc == expected_arg_num) {
		pars.par_file_name = argv[1];
		pars.seed = atof(argv[2]);
		return pars;
	} else {
		err("MISTAKE: there are 2 arguments: first - file name; second - seed ");
	};
	return pars;
};

void DoesInitiationPart(Parameters & pars) {
	string full_par_file_name = pars.par_file_name;

	if (FileExists(full_par_file_name.c_str())) {
		GetPars(pars);
	} else {
		cout << "file with the name "<<  pars.par_file_name <<
						" doesn't exist"<< endl;
	};

	pars.vars.gen.seed(pars.seed);

	InitLoop(pars);
};

//initializes parameters
void GetPars(Parameters & pars) {
  boost::property_tree::ptree pt;
  boost::property_tree::info_parser::read_info(pars.par_file_name, pt);
  // begin* Name a directory with etalon_mut_vec.dat (etalon mutation vector)
  pars.catalog_file_name = pt.get<string>
		("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file");
  pars.comparison_file_name = pt.get<string>
    ("Parameters.ProbeComparisonParameters.comparison_file_name");
  pars.number_of_parameters = 3;
	pars.algo.final_epsilon = pt.get<float>
	 ("algorithm_parameters.final_epsilon");
  pars.vars.max_dist = 1000;
  pars.algo.parameter_population_num = pt.get<float>
		("algorithm_parameters.parameter_population_number");
  pars.algo.parameter_population_num_alpha =
		pars.algo.parameter_population_num/2;
  pars.algo.p_acc_min = pt.get<float>
	 ("algorithm_parameters.p_acc_min");

  pars.search_space.driver_adv.min_ = pt.get<float>
		("search_space.driver_advantage.min");
  pars.search_space.driver_adv.max_ = pt.get<float>
		("search_space.driver_advantage.max");
  pars.search_space.driver_adv.true_val = pt.get<float>
		("search_space.driver_advantage.true");

  pars.search_space.mut_rate.min_ = pt.get<float>
		("search_space.mutation_rate.min");
  pars.search_space.mut_rate.max_ =  pt.get<float>
		("search_space.mutation_rate.max");
  pars.search_space.mut_rate.true_val = pt.get<float>
		("search_space.mutation_rate.true");

  pars.search_space.driver_mut_rate.min_ = pt.get<float>
		("search_space.driver_mutation_rate.min");
  pars.search_space.driver_mut_rate.max_ = pt.get<float>
		("search_space.driver_mutation_rate.max");
  pars.search_space.driver_mut_rate.true_val = pt.get<float>
    ("search_space.driver_mutation_rate.true");

	float stoping_level = 0.0001;
  float driver_adv_range = pars.search_space.driver_adv.max_-
                          pars.search_space.driver_adv.min_;
  float mut_rate_range = pars.search_space.mut_rate.max_-
                          pars.search_space.mut_rate.min_;
  float driver_mut_rate_range = pars.search_space.driver_mut_rate.max_-
                          pars.search_space.driver_mut_rate.min_;

  if (driver_adv_range/pars.search_space.driver_adv.min_ <= stoping_level) {
		pars.search_space.driver_adv.is_mutable = false;
	} else {
		pars.search_space.driver_adv.is_mutable = true;
	}
  if (mut_rate_range/pars.search_space.mut_rate.min_ <= stoping_level) {
		pars.search_space.mut_rate.is_mutable = false;
  } else {
    pars.search_space.mut_rate.is_mutable = true;
  };
  if (driver_mut_rate_range/pars.search_space.driver_mut_rate.min_ <= stoping_level) {
		pars.search_space.driver_mut_rate.is_mutable = false;
  } else {
  	pars.search_space.driver_mut_rate.is_mutable =true;
  };

  pars.vars.p_acc = 1000;
};


void MainLoop(Parameters & pars,
							std::vector<TreeComparison>& tree_comparison_vec)	{
#if defined(DEBUG_REGIME)
	std::cout<<"We are in main loop"<<std::endl;
#endif
	unsigned int iteration_num = 1;//is used for names of files
	pars.vars.max_dist = std::max_element(tree_comparison_vec.begin(),
	/* save current max_dist for inner function;*/tree_comparison_vec.end(),
	/* will be changed during exec*/comparison_of_two_trees)->results.tree_dist;
	SetTreeComparisonVec2File(tree_comparison_vec, pars, 0);
	while ((pars.vars.max_dist > pars.algo.final_epsilon) &&
				 (pars.vars.p_acc > pars.algo.p_acc_min)){//the begin of main loop
		std::sort(tree_comparison_vec.begin(), tree_comparison_vec.end(),
							comparison_of_two_trees);
	  WeightNormalization(tree_comparison_vec, 0,/*partial weight normalization */
												pars.algo.parameter_population_num_alpha);
		SetTreeComparisonVec2File(tree_comparison_vec, pars,
															iteration_num);
		unsigned int m_size = pars.number_of_parameters;//matrix_size

		ublas::vector<float> weighted_mean_vector;
		weighted_mean_vector = GetWeightedMeanVec(tree_comparison_vec, pars);

		ublas::matrix<float, ublas::row_major> epsilon_matrix(m_size, m_size);
		epsilon_matrix = ublas::zero_matrix<float>(m_size, m_size);
		epsilon_matrix = GetTwiceWeightedEmpiricalCovariance(
												weighted_mean_vector, tree_comparison_vec, pars);
#if defined(DEBUG_REGIME)
    cout<<"In func MainLoop"<<std::endl;
    cout<<"weighted_mean_vector="<<weighted_mean_vector<<endl;
    cout<<"epsilon_matrix="<<epsilon_matrix<<endl;
		int i;
    cin>>i;
#endif
		pars.vars.n_trials = 0;
		InnerLoop(pars, tree_comparison_vec, epsilon_matrix);
    float diff = tree_comparison_vec.size() -
                 pars.algo.parameter_population_num_alpha;
 		pars.vars.p_acc = (diff)/(pars.vars.n_trials*1.0);
		WeightNormalization(tree_comparison_vec,
												 pars.algo.parameter_population_num_alpha,
												 tree_comparison_vec.size());
		WeightNormalization(tree_comparison_vec,0,
												 tree_comparison_vec.size());

		iteration_num++;
		SetTreeComparisonVec2File(tree_comparison_vec, pars, iteration_num);
	};
};
