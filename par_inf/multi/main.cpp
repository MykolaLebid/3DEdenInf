//i/o includes
//#define __USE_MINGW_ANSI_STDIO 0
//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <fstream>

//parallel calculations
#include <omp.h>

//C++ 11
#include <random>
#include <algorithm> // std::sort
#include <vector>

//multivariate_gaussian_distribution
#include "multivariate_gaussian_distribution.h"
#include "cholesky.hpp"
#include "inverse_matrix.h"

//boost
//#include <boost/limits.hpp>
namespace ublas = boost::numeric::ublas;

//debugging
//#include "debug_functions.h"

//#define DEBUG_REGIME  // check regime
#include "prog_parameters.h"

void InitPars(const int argc, char *argv[], AllSimParameters & conf_settings); // function in main
void InitLoop(const AllSimParameters & conf_settings, Vars & var); // function in main
void InitWeight(const int parameter_population_num,
								std::vector <TreeComparison> & tree_comparison_vec); // function in main
//// next key function after main()
void MainLoop(const AllSimParameters & conf_settings, Vars & vars,
							std::vector<TreeComparison> & tree_comparison_vec);
////////////////////////////////////////////////////////////////////////////////
//// support functions for the ABC algorithm (in MainLoop())
ublas::vector<float> GetWeightedMeanVec(
	const std::vector<TreeComparison> & tree_comparison_vec,
	const AllSimParameters & conf_settings);
ublas::vector<float> op(const ublas::vector<float> & weighted_vector,
												TreeComparison T); // (in GetWeightedMeanVec())
TreeComparison GetSampleVarianceParDist(const std::vector<TreeComparison> &
																							tree_comparison_vec);
TreeComparison GetMeanParDist(const AllSimParameters & conf_settings,
															const std::vector<TreeComparison>& tree_comparison_vec);
void WeightNormalization(const int begin_index, const int end_index,
												 std::vector<TreeComparison>& tree_comparison_vec);
ublas::matrix<float, ublas::row_major> GetTwiceWeightedEmpiricalCovariance(
										const ublas::vector<float> & weighted_mean_vector,
										const std::vector<TreeComparison> & array_of_tree_comparisons,
										const AllSimParameters & conf_settings);
//// function in GetTwiceWeightedEmpiricalCovariance()
float SumOfSquareWeights(const std::vector<TreeComparison> & array_of_tree_comparisons,
												 const int parameter_population_num_alpha);
////
//////////////////////////////////////////////////////////////////////////////////
//
//// next key function in MainLoop()
void InnerLoop(const AllSimParameters & conf_settings, Vars & vars,
							 std::vector<TreeComparison> & tree_comparison_vec,
							 const ublas::matrix<float, ublas::row_major> & epsilon_matrix);
//////////////////////////////////////////////////////////////////////////////////
//// support functions in InnerLoop()
void SetNewWeights2NewTreeComparisons(
												const ublas::matrix<float, ublas::row_major> & epsilon_matrix,
												const AllSimParameters & conf_settings,
												std::vector<TreeComparison>& tree_comparison_vec);
//// support functions in SetNewWeights2NewTreeComparisons())
float GetGaussianKernel(
		const std::vector<TreeComparison> & array_of_tree_comparisons,
		const ublas::matrix<float, ublas::row_major> & epsilon_matrix,
		const int i, const int j, const AllSimParameters & conf_settings);
////
//////////////////////////////////////////////////////////////////////////////////
//
//// next key function in InnerLoop()
void SubInnerLoop(const std::vector<TreeComparison> & tree_comparison_vec,
					        const AllSimParameters & conf_settings, Vars & var,
					        ublas::matrix<float, ublas::row_major> & Cholesky_matrix);
//////////////////////////////////////////////////////////////////////////////////
//// support functions in SubInnerLoop()
int GetIndexOfParameterVector(
		const std::vector<TreeComparison> array_of_tree_comparisons,
		const int parameter_population_num_alpha,
		std::mt19937 & gen);
ublas::vector<float> GetMultiNormalDistSample(
                        const ublas::vector<float> & parameter_vector,
                        const ublas::matrix<float, ublas::row_major> & L,
                        const AllSimParameters & conf_settings, Vars & var);
// support functions in GetIndexOfParameterVector()
ublas::vector<float> ReflectedCoords(const AllSimParameters & conf_settings,
																		 ublas::vector<float> & new_parameter_vector);
//// support functions in GetReflectedCoord()
float GetReflectedCoord(const float old_coord, const float l_b/*left border*/,
																						   const float r_b/*left border*/);
//////// support functions in GetReflectedCoord
void SetProjs(const std::vector<TreeComparison> & array_of_tree_comparisons,
							const std::vector<bool> is_par_vec,const int i, const int j,
							ublas::vector<float> & proj_i_vector,
							ublas::vector<float> & proj_j_vector);
////
//////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]){
  AllSimParameters conf_settings;
  Vars vars;
  // the number of arguments in command line
  const int expected_arg_num = 4;
  if (argc == expected_arg_num){
		InitPars(argc, argv, conf_settings);
	} else {
		err("MISTAKE: there are 3 arguments: first - related path to a configuration (parameter) file; second - name of configuration (parameter) file; third  - seed ");
	};

	std::string full_par_file_name =
						conf_settings.basic_sim_parameters.related_path_to_par_file +
						conf_settings.basic_sim_parameters.par_file_name;
	if (FileExists(full_par_file_name.c_str())){
		vars.gen.seed(conf_settings.basic_sim_parameters.seed);
		InitLoop(conf_settings, vars);
		std::vector<TreeComparison> tree_comparison_vec;
		int par_pop_num = conf_settings.inference_parameters.
												algorithm_parameters.parameter_population_number;
		GetTreeComparisonVecFromFile(conf_settings, par_pop_num, tree_comparison_vec);
    InitWeight(par_pop_num, tree_comparison_vec);
#if defined(DEBUG_REGIME)
		SetTreeComparisonVec2File(tree_comparison_vec, conf_settings, 0);
#endif
		MainLoop(conf_settings, vars, tree_comparison_vec);
	} else {
		std::cout << "file with the name "<<  full_par_file_name <<
						" doesn't exist"<< std::endl;
	};
  return 0;
};

void SetProjs(const std::vector<TreeComparison> & array_of_tree_comparisons,
							const std::vector <bool> is_par_vec, const int i, const int j,
							ublas::vector<float> & proj_i_vector,
							ublas::vector<float> & proj_j_vector) {
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


float GetReflectedCoord(const float old_coord, const float l_b,const float r_b) {
	//std::cout<<"old_coord="<<old_coord<<"left_border="<<l_b<<"right_border="<<r_b<<endl;
	float diff_min = old_coord - l_b;
	if (diff_min<0) {
		return (l_b + abs(diff_min));
	} else {
		float diff_max = old_coord - r_b;
		if (diff_max<0) err("GetReflectedCoord: diff_max<0 ");
		return  (r_b - diff_max);
	};
};

// reflect coords of parameters
// until all of them will be in acceptable parameter range
ublas::vector<float> ReflectedCoords(const AllSimParameters & conf_settings,
																		 ublas::vector<float> & new_parameter_vector){
  std::vector <bool> is_in_dom_vec = {false, false, false};
  do {
		if ((GetLeftBorderDriverAdvantage(conf_settings)<=new_parameter_vector[0] &&
			GetRigthBorderDriverAdvantage(conf_settings)>=new_parameter_vector[0]) || is_in_dom_vec[0] ) {
			  is_in_dom_vec[0] = true;
		} else {
			//std::cout<<"ReflectedCoords driver_adv"<<std::endl;
			new_parameter_vector[0] = GetReflectedCoord(new_parameter_vector[0],
															   GetLeftBorderDriverAdvantage(conf_settings),
															   GetRigthBorderDriverAdvantage(conf_settings));
		}; // process driver_adv coord
		if ((GetLeftBorderMutationRate(conf_settings)<=new_parameter_vector[1] &&
				GetRigthBorderMutationRate(conf_settings)>=new_parameter_vector[1]) || is_in_dom_vec[1]) {
			is_in_dom_vec[1] = true;
		} else {
			//std::cout<<"ReflectedCoords mut_rate"<<std::endl;
			new_parameter_vector[1] = GetReflectedCoord(new_parameter_vector[1],
															   GetLeftBorderMutationRate(conf_settings),
															   GetRigthBorderMutationRate(conf_settings));
		}; // process mut_rate coord
		if ((GetLeftBorderDriverMutationRate(conf_settings)<=new_parameter_vector[2] &&
				GetRigthBorderDriverMutationRate(conf_settings)>=new_parameter_vector[2]) || is_in_dom_vec[2]) {
			is_in_dom_vec[2] = true;
		} else {
			//std::cout<<"ReflectedCoords driver_mut_rate"<<std::endl;
			new_parameter_vector[2] = GetReflectedCoord(new_parameter_vector[2],
															   GetLeftBorderDriverMutationRate(conf_settings),
															   GetRigthBorderDriverMutationRate(conf_settings));
		}; // process driver_mut_rate coord
	} while (is_in_dom_vec[0] == false ||
					 is_in_dom_vec[1] == false ||
					 is_in_dom_vec[2] == false );
	return new_parameter_vector;
};
//    return parameter_vector;


ublas::vector<float> GetMultiNormalDistSample(
                        const ublas::vector<float> & parameter_vector,
                        const ublas::matrix<float, ublas::row_major> & L,
                        const AllSimParameters & conf_settings, Vars & vars) {
	int size_matrix = GetSearchParNumber();
	ublas::vector<float> new_parameter_vector(size_matrix);
	std::normal_distribution<float> n_d(0,1);
	for (int i=0; i < size_matrix; i++) new_parameter_vector[i] = n_d(vars.gen);
	new_parameter_vector = prod(new_parameter_vector, L);
#if defined(DEBUG_REGIME)
	std::cout<<"We are in GetMultiNormalDistSample"<<std::endl;
  std::cout<<"before adding present param new_parameter_vector="<<
		new_parameter_vector<<std::endl;
std::cout<<"present param: ="<<
		parameter_vector<<std::endl;
#endif

	if (conf_settings.inference_parameters.
				log_scale_search_space.driver_advantage == 0) new_parameter_vector[0] = 0;
	if (conf_settings.inference_parameters.
				log_scale_search_space.mutation_rate == 0) new_parameter_vector[1] = 0;
	if (conf_settings.inference_parameters.
				log_scale_search_space.driver_mutation_rate == 0) new_parameter_vector[2] = 0;

	new_parameter_vector += parameter_vector;

	// reflect coords of parameters
	// until all of them will be in acceptable parameter range
	ReflectedCoords(conf_settings, new_parameter_vector);
#if defined(DEBUG_REGIME)
	std::cout<<"After reflection new_parameter_vector="<<
		new_parameter_vector<<std::endl;
#endif

	return new_parameter_vector;
};

int GetIndexOfParameterVector (
		const std::vector<TreeComparison> array_of_tree_comparisons,
		const int parameter_population_num_alpha,
		std::mt19937 & gen) {
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

float GetGaussianKernel(
				const std::vector<TreeComparison> & array_of_tree_comparisons,
				const ublas::matrix<float, ublas::row_major> & epsilon_matrix,
				const int i, const int j, const AllSimParameters & conf_settings) {
#if defined(DEBUG_REGIME)
  std::cout<<"In GetGaussianKernel="<<std::endl;
#endif
//	int parameter_number;
	std::vector <bool> is_par_vec;
	if (conf_settings.inference_parameters.log_scale_search_space.driver_advantage == 0) {
		is_par_vec.push_back(false);
	} else {
		is_par_vec.push_back(true);
	};
	if (conf_settings.inference_parameters.log_scale_search_space.mutation_rate == 0) {
		is_par_vec.push_back(false);
	} else {
		is_par_vec.push_back(true);
	};
	if (conf_settings.inference_parameters.log_scale_search_space.driver_mutation_rate == 0) {
		is_par_vec.push_back(false);
	} else {
		is_par_vec.push_back(true);
	};

	int num_par = 0;
	for (unsigned int it = 0; it < is_par_vec.size(); it++) if (is_par_vec[it]) num_par++;
	range r(0, num_par);
#if defined(DEBUG_REGIME)
	std::cout << "num_par="<<num_par <<std::endl;
#endif
	ublas::matrix<float, ublas::row_major> proj_epsilon_matrix(num_par, num_par),
																	inverse_proj_epsilon_matrix(num_par, num_par);
	ublas::vector<float> proj_i_vector(num_par), proj_j_vector(num_par),
					             diff_proj_i_j_vector(num_par), result_of_mult(num_par);
  SetProjs(array_of_tree_comparisons, is_par_vec, i, j, proj_i_vector, proj_j_vector);

	diff_proj_i_j_vector = proj_i_vector - proj_j_vector;
#if defined(DEBUG_REGIME)
  std::cout<<"proj_i_vector="<<proj_i_vector<<std::endl;
  std::cout<<"proj_j_vector="<<proj_j_vector<<std::endl;
  std::cout<<"diff_proj_i_j_vector="<<diff_proj_i_j_vector<<std::endl;
#endif

  int index_row = 0;
  for(int row = 0; row < GetSearchParNumber(); row++) {
		if (is_par_vec[row]) {
			int index_col = 0;
			for(int col = 0; col < GetSearchParNumber(); col++) {
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
	//int i_00;
	//cin>>i_00;
#endif
	return inner_prod(diff_proj_i_j_vector, result_of_mult);
};

void SubInnerLoop(const std::vector<TreeComparison> & tree_comparison_vec,
									const AllSimParameters & conf_settings, Vars & vars,
									ublas::matrix<float, ublas::row_major> & Cholesky_matrix){
	int matrix_size = GetSearchParNumber();
	int parameter_population_num_alpha = conf_settings.inference_parameters.
														algorithm_parameters.parameter_population_num_alpha;
	int index = GetIndexOfParameterVector(tree_comparison_vec,
																			parameter_population_num_alpha,
																			vars.gen);
	ublas::vector<float> parameter_vector(matrix_size);
	parameter_vector[0] = tree_comparison_vec[index].par_point.driver_adv;
	parameter_vector[1] = tree_comparison_vec[index].par_point.mut_rate;
	parameter_vector[2] = tree_comparison_vec[index].par_point.driver_mut_rate;

	parameter_vector = GetMultiNormalDistSample(parameter_vector, Cholesky_matrix,
																							conf_settings, vars);
  ParPoint par_point={parameter_vector [0],
											parameter_vector [1],
											parameter_vector [2]};
	SimulateAndWrite2File(par_point, conf_settings, GetRandNumber(vars));
	vars.n_trials++;
};

float SumOfSquareWeights(const std::vector<TreeComparison> & array_of_tree_comparisons,
												 const int parameter_population_num_alpha) {
  float sum = 0;
  for (int i = 0; i<parameter_population_num_alpha; i++)
		sum += array_of_tree_comparisons[i].results.weight*
					   array_of_tree_comparisons[i].results.weight;
  return sum;
};

void SetNewWeights2NewTreeComparisons(
												const ublas::matrix<float, ublas::row_major> & epsilon_matrix,
												const AllSimParameters & conf_settings,
												std::vector<TreeComparison>& tree_comparison_vec){
	int index_of_begin = conf_settings.inference_parameters.
											 algorithm_parameters.parameter_population_num_alpha;
	int index_of_end = tree_comparison_vec.size();
	for (int i = index_of_begin; i < index_of_end; i++) {
		float sum_denomination = 0;
		for(int j = 0; j < index_of_begin; j++) {
			sum_denomination += tree_comparison_vec[j].results.weight *
			GetGaussianKernel(tree_comparison_vec, epsilon_matrix, i, j, conf_settings);
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

void InnerLoop(const AllSimParameters & conf_settings, Vars & vars,
							 std::vector<TreeComparison> & tree_comparison_vec,
							 const ublas::matrix<float, ublas::row_major> & epsilon_matrix){
	int matrix_size = GetSearchParNumber();
	int parameter_population_number = conf_settings.inference_parameters.
														algorithm_parameters.parameter_population_number;
	int parameter_population_num_alpha = conf_settings.inference_parameters.
														algorithm_parameters.parameter_population_num_alpha;
	vars.max_dist = tree_comparison_vec[
								tree_comparison_vec.size() - parameter_population_num_alpha - 1
																			].results.tree_dist;
  // init file for new loop
  std::string file_name(GetFullComparisonFileName(conf_settings));
	remove(file_name.c_str());  // end init file for new loop
  //Cholesky_decompose
	ublas::matrix<float, ublas::row_major> L (matrix_size, matrix_size);
	L = ublas::zero_matrix<float>(matrix_size, matrix_size);
	cholesky_decompose(epsilon_matrix, L);
	// end Cholesky_decompose
#if defined(DEBUG_REGIME)
	std::cout<<"In InnerLoop"<<std::endl;
	std::cout<<"cholesky_decompose="<<L<<std::endl;
#endif

  int diff = parameter_population_number - parameter_population_num_alpha;
	do{//shared(flag) #pragma omp parallel for
			int boundary = 2.0 * parameter_population_number;
  		//download new mingw in a case of problems for OpenMP
      #pragma omp parallel for
			for(int i = 0; i < boundary; i++) {
				SubInnerLoop(tree_comparison_vec, conf_settings, vars, L);
			};
	} while(GetNumLinesInComparisonFile(conf_settings) <= diff);

  std::vector<TreeComparison> tree_comparison_vec_plus;

  GetTreeComparisonVecFromFile(conf_settings, diff, tree_comparison_vec_plus);
	//send results to main array array_of_tree_comparisons from additional array
	for(int i = parameter_population_num_alpha;
			i < parameter_population_number;
			i++) tree_comparison_vec[i] = tree_comparison_vec_plus[i -
		                           parameter_population_num_alpha];
	SetNewWeights2NewTreeComparisons(epsilon_matrix, conf_settings, tree_comparison_vec);
};

ublas::matrix<float, ublas::row_major> GetTwiceWeightedEmpiricalCovariance(
										const ublas::vector<float> & weighted_mean_vector,
										const std::vector<TreeComparison> & array_of_tree_comparisons,
										const AllSimParameters & conf_settings) {
  int size_matrix = GetSearchParNumber();
	int parameter_population_num_alpha = conf_settings.inference_parameters.
													 algorithm_parameters.parameter_population_num_alpha;
  //init epsion_matrix
  ublas::matrix<float, ublas::row_major> epsilon_matrix(size_matrix, size_matrix);
  epsilon_matrix = ublas::zero_matrix<float>(size_matrix, size_matrix);
#if defined(DEBUG_REGIME)
    std::cout<<"In function GetTwiceWeightedEmpiricalCovariance:"<<std::endl;
#endif

  for (int i = 0; i < parameter_population_num_alpha; i++)	{
    ublas::vector<float> parameter_vector(size_matrix);

    parameter_vector(0) = array_of_tree_comparisons[i].par_point.driver_adv;
    parameter_vector(1) = array_of_tree_comparisons[i].par_point.mut_rate;
    parameter_vector(2) = array_of_tree_comparisons[i].par_point.driver_mut_rate;

    ublas::vector<float> diff_mean_vector(size_matrix);

    diff_mean_vector(0) = parameter_vector[0] - weighted_mean_vector[0];
    diff_mean_vector(1) = parameter_vector[1] - weighted_mean_vector[1];
    diff_mean_vector(2) = parameter_vector[2] - weighted_mean_vector[2];
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

 // for Choletski decomposition (create positive definite matrices)
 // will be neglected in
  for (unsigned int i = 0;i < epsilon_matrix.size1();i++)
		if (epsilon_matrix(i,i)==0) epsilon_matrix(i,i) = 1e-12;

  for(unsigned int i = 0; i < epsilon_matrix.size1(); i++)
		for(unsigned int j = 0; j < epsilon_matrix.size2(); j++)
				if (epsilon_matrix(i,j)==0) epsilon_matrix(i,j) = 1e-15;

  return epsilon_matrix;
};

void WeightNormalization(const int begin_index, const int end_index,
												 std::vector<TreeComparison> & tree_comparison_vec){
	float weight_sum_alpha = 0;
	for(int i = begin_index; i < end_index; i++)
		weight_sum_alpha += tree_comparison_vec[i].results.weight;
 // cout<< "parameter_population_num_alpha="<<parameter_population_num_alpha;
 // cout<< "weight_sum_alpha="<<weight_sum_alpha;
	for(int i = begin_index; i < end_index; i++)
		tree_comparison_vec[i].results.weight =
			((tree_comparison_vec[i].results.weight * 1.0) / (1.0 * weight_sum_alpha));
};

TreeComparison GetMeanParDist(const AllSimParameters & conf_settings,
															const std::vector<TreeComparison>& tree_comparison_vec){
	TreeComparison A; // 1) sum up all statistical results;
	                  // 2) find the average;
	                  // 3) find relative error with true results;
  // 1)
	A.par_point.driver_adv = 0;
	A.par_point.mut_rate = 0;
	A.par_point.driver_mut_rate = 0;
	unsigned int vec_size = tree_comparison_vec.size();
	for (unsigned int i = 0; i < vec_size; i++) {
		A.par_point.driver_adv += tree_comparison_vec[i].par_point.driver_adv;
		A.par_point.mut_rate += tree_comparison_vec[i].par_point.mut_rate;
		A.par_point.driver_mut_rate += tree_comparison_vec[i].par_point.driver_mut_rate;
	};
	// 2)
	A.par_point.driver_adv = A.par_point.driver_adv/(1.0*vec_size);
	A.par_point.mut_rate = A.par_point.mut_rate/(1.0*vec_size);
	A.par_point.driver_mut_rate = A.par_point.driver_mut_rate/(1.0*vec_size);
	// 3)
	float true_driver_adv = conf_settings.basic_sim_parameters.evolution.driver_adv;
	A.par_point.driver_adv = (true_driver_adv - A.par_point.driver_adv) / true_driver_adv;
	float true_mut_rate = conf_settings.basic_sim_parameters.evolution.mutation_rate;
	A.par_point.mut_rate = (true_mut_rate - A.par_point.mut_rate) / true_mut_rate;
	float true_dr_mut_rate = conf_settings.basic_sim_parameters.evolution.driver_mutation_rate;
	A.par_point.driver_mut_rate = (true_dr_mut_rate - A.par_point.driver_mut_rate)/
	                              true_dr_mut_rate;
	return A;
};

TreeComparison GetSampleVarianceParDist(const std::vector<TreeComparison> &
																							tree_comparison_vec){
	TreeComparison A; // 1) sum up all statistical results;
	                  // 2) find the average;
	                  // 3) sample variation;
  // 1)
	A.par_point.driver_adv = 0;
	A.par_point.driver_mut_rate = 0;
	A.par_point.mut_rate = 0;
	unsigned int vec_size = tree_comparison_vec.size();
	for (unsigned int i = 0; i < vec_size; i++) {
		A.par_point.driver_adv += tree_comparison_vec[i].par_point.driver_adv;
		A.par_point.mut_rate += tree_comparison_vec[i].par_point.mut_rate;
		A.par_point.driver_mut_rate += tree_comparison_vec[i].par_point.driver_mut_rate;
	};

	// 2)
	A.par_point.driver_adv = A.par_point.driver_adv/(1.0*vec_size);
	A.par_point.mut_rate = A.par_point.mut_rate/(1.0*vec_size);
	A.par_point.driver_mut_rate = A.par_point.driver_mut_rate/(1.0*vec_size);

	// 3)
  float var_1 = 0, var_2 = 0, var_3 = 0;
	for (unsigned int i = 0; i < vec_size; i++) {
		var_1 += (A.par_point.driver_adv - tree_comparison_vec[i].par_point.driver_adv) *
						 (A.par_point.driver_adv - tree_comparison_vec[i].par_point.driver_adv);
		var_2 += (A.par_point.mut_rate - tree_comparison_vec[i].par_point.mut_rate)*
		         (A.par_point.mut_rate - tree_comparison_vec[i].par_point.mut_rate);
		var_3 += (A.par_point.driver_mut_rate - tree_comparison_vec[i].par_point.driver_mut_rate) *
						 (A.par_point.driver_mut_rate - tree_comparison_vec[i].par_point.driver_mut_rate);
	};

	A.par_point.driver_adv = var_1/(1.0 * vec_size);
	A.par_point.mut_rate = var_2/(1.0 * vec_size);
	A.par_point.driver_mut_rate = var_3/(1.0 * vec_size);

	return A;
};

ublas::vector<float> op(const ublas::vector<float> & weighted_vector,
												TreeComparison T) {
    ublas::vector<float> w_v(weighted_vector.size());
    w_v(0) = weighted_vector(0) + T.results.weight * T.par_point.driver_adv;
    w_v(1) = weighted_vector(1) + T.results.weight * T.par_point.mut_rate;
    w_v(2) = weighted_vector(2) + T.results.weight * T.par_point.driver_mut_rate;
    return w_v;
};

ublas::vector<float> GetWeightedMeanVec(
	const std::vector<TreeComparison> & tree_comparison_vec,
	const AllSimParameters & conf_settings){

  int size_of_matrix = GetSearchParNumber();
  int parameter_population_num_alpha = conf_settings.inference_parameters.
													algorithm_parameters.parameter_population_num_alpha;
  ublas::vector<float> weighted_mean_vector(size_of_matrix);

  weighted_mean_vector = ublas::zero_vector<float> (size_of_matrix);
  weighted_mean_vector = std::accumulate(
                            tree_comparison_vec.begin(),
                            tree_comparison_vec.begin() +
																				parameter_population_num_alpha,
                            weighted_mean_vector,
					                  op);
	if (conf_settings.inference_parameters.log_scale_search_space.driver_advantage == 0) {
		weighted_mean_vector[0] = conf_settings.basic_sim_parameters.evolution.driver_adv;
	};

	if (conf_settings.inference_parameters.log_scale_search_space.mutation_rate == 0) {
		weighted_mean_vector[1] = conf_settings.basic_sim_parameters.evolution.mutation_rate;
	};

	if (conf_settings.inference_parameters.log_scale_search_space.driver_mutation_rate == 0) {
		weighted_mean_vector[2] = conf_settings.basic_sim_parameters.evolution.driver_mutation_rate;
	};

  return weighted_mean_vector;
};
//
void MainLoop(const AllSimParameters & conf_settings, Vars & vars,
							std::vector<TreeComparison> & tree_comparison_vec){
#if defined(DEBUG_REGIME)
	std::cout<<"We are in the main loop"<<std::endl;
#endif
	unsigned int iteration_num = 1; //is used for names of files
	vars.max_dist = std::max_element(tree_comparison_vec.begin(),
	/* save current max_dist for inner function;*/ tree_comparison_vec.end(),
	/* will be changed during exec*/CompareTwoTrees)->results.tree_dist;
#if defined(DEBUG_REGIME)
	SetTreeComparisonVec2File(tree_comparison_vec, conf_settings, 0);
#endif
	float final_epsilon = conf_settings.inference_parameters.
													algorithm_parameters.final_epsilon;
	float p_acc_min = conf_settings.inference_parameters.
													algorithm_parameters.p_acc_min;
	vars.p_acc = 1000;
	float final_iteration = conf_settings.inference_parameters.
													algorithm_parameters.final_iteration;
	int parameter_population_num_alpha = conf_settings.inference_parameters.
														algorithm_parameters.parameter_population_num_alpha;

	while ((vars.max_dist > final_epsilon) && (vars.p_acc > p_acc_min) &&
				 (iteration_num < final_iteration)){//the begin of main loop
		std::sort(tree_comparison_vec.begin(), tree_comparison_vec.end(),
							CompareTwoTrees);
	  WeightNormalization(0, parameter_population_num_alpha, tree_comparison_vec);
#if defined(DEBUG_REGIME)
		SetTreeComparisonVec2File(tree_comparison_vec, conf_settings, iteration_num);
#endif
		unsigned int m_size = GetSearchParNumber();//matrix_size

		ublas::vector<float> weighted_mean_vector;

		weighted_mean_vector = GetWeightedMeanVec(tree_comparison_vec, conf_settings);

		ublas::matrix<float, ublas::row_major> epsilon_matrix(m_size, m_size);
		epsilon_matrix = ublas::zero_matrix<float>(m_size, m_size);
		epsilon_matrix = GetTwiceWeightedEmpiricalCovariance(
												weighted_mean_vector, tree_comparison_vec, conf_settings);
#if defined(DEBUG_REGIME)
    std::cout<< "In func MainLoop" <<std::endl;
    std::cout<< "weighted_mean_vector=" <<weighted_mean_vector<<std::endl;
    std::cout<< "epsilon_matrix=" <<epsilon_matrix<<std::endl;
		//int i;
    //std::cin>>i;
#endif
		vars.n_trials = 0;
		InnerLoop(conf_settings, vars, tree_comparison_vec, epsilon_matrix);
    float diff = tree_comparison_vec.size() - parameter_population_num_alpha;
 		vars.p_acc = (diff)/(vars.n_trials * 1.0);
		WeightNormalization(parameter_population_num_alpha,
												 tree_comparison_vec.size(), tree_comparison_vec);
		WeightNormalization(0, tree_comparison_vec.size(), tree_comparison_vec);
		iteration_num ++;
		TreeComparison treeComparisonRelativeError =
																			GetMeanParDist(conf_settings, tree_comparison_vec);
		TreeComparison treeComparisonSampleVariance =
														GetSampleVarianceParDist(tree_comparison_vec);
		AddToLogFile(conf_settings, iteration_num, vars.n_trials,
								 treeComparisonRelativeError,
								 treeComparisonSampleVariance);
#if defined(DEBUG_REGIME)
		SetTreeComparisonVec2File(tree_comparison_vec, conf_settings, iteration_num);
#endif

	};
};

void InitWeight(const int parameter_population_num,
								std::vector <TreeComparison> & tree_comparison_vec){
	for(int i=0; i < parameter_population_num; i++) {
		tree_comparison_vec[i].results.weight = 1.0/parameter_population_num;
	};
};

void InitLoop(const AllSimParameters & conf_settings, Vars & vars){
	std::string full_comparison_file_name = GetFullComparisonFileName(conf_settings);
	remove(full_comparison_file_name.c_str());
	std::string full_log_file_name = GetPathToEtalonResultCatalogue(conf_settings) +
					conf_settings.inference_parameters.basic_settings.log_file_name;
	remove(full_log_file_name.c_str());
#if defined(DEBUG_REGIME)
	std::cout <<"Search space:"<<"["<<GetLeftBorderDriverAdvantage(conf_settings)<<","<<GetRigthBorderDriverAdvantage(conf_settings)<< "]"<<std::endl<<
						  "["<<GetLeftBorderMutationRate(conf_settings)<<","<<GetRigthBorderMutationRate(conf_settings)<< "]"<<std::endl<<
						  "["<<GetLeftBorderDriverMutationRate(conf_settings)<<","<<GetRigthBorderDriverMutationRate(conf_settings)<< "]"<<std::endl;
#endif // defined

	std::uniform_real_distribution<float> driver_adv_gen(GetLeftBorderDriverAdvantage(conf_settings),
																											 GetRigthBorderDriverAdvantage(conf_settings));
	std::uniform_real_distribution<float> mut_rate_gen(GetLeftBorderMutationRate(conf_settings),
																										 GetRigthBorderMutationRate(conf_settings));
	std::uniform_real_distribution<float> driver_mut_rate_gen(GetLeftBorderDriverMutationRate(conf_settings),
																														GetRigthBorderDriverMutationRate(conf_settings));
	int parameter_population_number = conf_settings.inference_parameters.
																		  algorithm_parameters.parameter_population_number;
	do {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  		  #pragma omp parallel for
			  for(int i = 0; i <= parameter_population_number; i++) {
					SimulateAndWrite2File({driver_adv_gen(vars.gen),
																 mut_rate_gen(vars.gen),
																 driver_mut_rate_gen(vars.gen)},
																 conf_settings, GetRandNumber(vars));
			  };
	} while (GetNumLinesInComparisonFile(conf_settings) <= parameter_population_number);
};

void InitPars(const int argc, char *argv[], AllSimParameters & conf_settings) {
		conf_settings.basic_sim_parameters.related_path_to_par_file = argv[1];
		conf_settings.basic_sim_parameters.par_file_name = argv[2];
		GetAllParameters(conf_settings.basic_sim_parameters,
										 conf_settings.inference_parameters);
		conf_settings.basic_sim_parameters.seed = atoi(argv[3]);
};
