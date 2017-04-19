
//i/o includes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//#include <boost/filesystem.hpp>

//parallel calculations
#include <omp.h>

//
#include <random>
#include <algorithm>    // std::sort
#include <vector>
#include <ctime>
//#include "stdafx.h"
#include <fstream>

#include <boost/limits.hpp>
// /multivariate_gaussian_distribution/ //
#include "multivariate_gaussian_distribution.h"
#include "cholesky.hpp"
#include "inverse_matrix.h"

#include "debug_functions.h"
#define CHECK_REGIME  // check regime





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

bool comparison_of_two_trees(TreeComparison tree_results_1,
														 TreeComparison tree_results_2) {
    return (tree_results_1.results.tree_dist <
						  tree_results_2.results.tree_dist);
};

// the file is named by "%s/file_tree_comparison.dat" from cancer_3d
void SimulateAndWrite2File(const ParPoint par_point, Parameters & pars) {
    char command_line[256];
    int rand = (par_point.driver_adv +
								par_point.driver_mut_rate +
								par_point.mut_rate) * 1000000000;
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
    cout <<"command_line[256]="<<command_line<<endl;
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
													 >> A.par_point.driver_mut_rate
													 >> A.par_point.mut_rate
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
													 >> A.par_point.driver_mut_rate
													 >> A.par_point.mut_rate
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

	sprintf(name_file,"%s/%s_%d.dat", pars.catalog_file_name.c_str(),
					                          pars.comparison_file_name.c_str(),
					                          number_of_iteration);
	file_tree_comparison.open(name_file);

	if (file_tree_comparison.is_open()==false) cout<<"err..."<<endl;

	file_tree_comparison.precision(5);
	file_tree_comparison.setf(std::ios::fixed, std:: ios::floatfield);


	for(unsigned int i=0; i < tree_comparison_vec.size(); i++)
		file_tree_comparison<<"\t"<<
			tree_comparison_vec[i].par_point.driver_adv<<"\t"<<
			tree_comparison_vec[i].par_point.driver_mut_rate<<"\t"<<
			tree_comparison_vec[i].par_point.mut_rate<<"\t"<<
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
	std::uniform_real_distribution<float> driver_adv_gen(/*(0.001, 0.999)*/
																				pars.search_space.driver_adv.min_,
																				pars.search_space.driver_adv.max_ );
	std::uniform_real_distribution<float> mut_rate_gen(/*(0.001, 0.999)*/
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
	float weight_sum_alpha=0;
	for(int i = begin_index; i<end_index; i++)
		weight_sum_alpha += tree_comparison_vec[i].results.weight;

 // cout<< "parameter_population_num_alpha="<<parameter_population_num_alpha;
 // cout<< "weight_sum_alpha="<<weight_sum_alpha;
	for(int i=begin_index; i<end_index; i++)
		tree_comparison_vec[i].results.weight=
			((tree_comparison_vec[i].results.weight*1.0)/
				(1.0*weight_sum_alpha));
};

ublas::vector<float> op(const ublas::vector<float> & weighted_vector,
												TreeComparison T) {
    ublas::vector<float> w_v(weighted_vector.size());
    w_v(0) = weighted_vector(0) + T.results.weight*T.par_point.driver_adv;
    w_v(1) = weighted_vector(1) + T.results.weight*T.par_point.mut_rate;
    w_v(2) = weighted_vector(2) + T.results.weight*T.par_point.driver_mut_rate;
    return w_v;
};

float sum_of_square_weights(std::vector<TreeComparison> & array_of_tree_comparisons,
														int parameter_population_num_alpha) {
  float sum = 0;
  for (int i=0; i<parameter_population_num_alpha; i++)
		sum += array_of_tree_comparisons[i].results.weight*
					   array_of_tree_comparisons[i].results.weight;
  return sum;
};

ublas::vector<float> GetWeightedMeanVec(std::vector<TreeComparison> &
																					tree_comparison_vec,
																				int size_of_matrix,
																				int parameter_population_num_alpha) {
  ublas::vector<float> weighted_mean_vector(size_of_matrix);

  weighted_mean_vector = ublas::zero_vector<float> (size_of_matrix);
  weighted_mean_vector = std::accumulate(
                            tree_comparison_vec.begin(),
                            tree_comparison_vec.begin() +
																				parameter_population_num_alpha,
                            weighted_mean_vector,
                            op);
  return weighted_mean_vector;
};


ublas::matrix<float, ublas::row_major> get_twice_weighted_empirical_covariance(
												ublas::vector<float> & weighted_mean_vector,
												std::vector<TreeComparison>& array_of_tree_comparisons,
												int size_matrix, int parameter_population_num_alpha) {
  //init epsion_matrix
  ublas::matrix<float, ublas::row_major> epsilon_matrix(size_matrix, size_matrix);
  epsilon_matrix = ublas::zero_matrix<float>(size_matrix, size_matrix);

  for (int i=0; i<parameter_population_num_alpha; i++)	{
    ublas::vector<float> parameter_vector(size_matrix);

    parameter_vector(0) = array_of_tree_comparisons[i].par_point.driver_adv;
    parameter_vector(1) = array_of_tree_comparisons[i].par_point.mut_rate;
    parameter_vector(2) = array_of_tree_comparisons[i].par_point.driver_mut_rate;

    ublas::vector<float> diff_mean_vector(size_matrix);
    diff_mean_vector[0] = parameter_vector[0] - weighted_mean_vector[0];
    diff_mean_vector[1]=parameter_vector[1]-weighted_mean_vector[1];
    diff_mean_vector[2]=parameter_vector[2]-weighted_mean_vector[2];

    cout<<ublas::outer_prod(diff_mean_vector,diff_mean_vector)<<endl;
    epsilon_matrix += array_of_tree_comparisons[i].results.weight*
                       (ublas::outer_prod(diff_mean_vector,diff_mean_vector));
  };

  float sum_sq_wei = sum_of_square_weights(array_of_tree_comparisons,
																					 parameter_population_num_alpha);
  epsilon_matrix*=(2.0/(1-sum_sq_wei));
  cout<<"weighted_mean_vector="<<weighted_mean_vector<<endl;
  cout<<"epsilon_matrix="<<epsilon_matrix<<endl;
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


ublas::vector<float> GetMultiNormalDistSample(
												ublas::matrix<float, ublas::row_major> & L,
                        ublas::vector<float> & parameter_vector,
                        mt19937 & gen) {
	int size_matrix = L.size1();

	ublas::vector<float> new_parameter_vector(size_matrix);
	std::normal_distribution<float> n_d(0,1);

	new_parameter_vector[0] = n_d(gen);
	new_parameter_vector[1] = n_d(gen);
	new_parameter_vector[2] = n_d(gen);

	new_parameter_vector = prod( new_parameter_vector, L);
	new_parameter_vector += parameter_vector;

	new_parameter_vector[0] = std::abs(new_parameter_vector[0]);
	new_parameter_vector[1] = std::abs(new_parameter_vector[1]);
	new_parameter_vector[2] = std::abs(new_parameter_vector[2]);
	return new_parameter_vector;
};

float get_gaussian_kernel(
				std::vector<TreeComparison>& array_of_tree_comparisons,
				ublas::matrix<float, ublas::row_major> & epsilon_matrix,
				int i, int j) {
	int number_of_parameters = 3;
	range r(0, number_of_parameters);
	ublas::matrix<float,ublas::row_major>
		proj_epsilon_matrix(number_of_parameters, number_of_parameters),
		inverse_proj_epsilon_matrix(number_of_parameters, number_of_parameters);

	ublas::vector<float> proj_i_vector(number_of_parameters),
											 proj_j_vector(number_of_parameters),
					             diff_proj_i_j_vector(number_of_parameters),
					             result_of_mult(number_of_parameters);

	proj_i_vector(0) = array_of_tree_comparisons[i].par_point.driver_adv;
	proj_j_vector(0) = array_of_tree_comparisons[j].par_point.driver_adv;
	diff_proj_i_j_vector = proj_i_vector - proj_j_vector;

	proj_epsilon_matrix = project(epsilon_matrix, r, r);
	InvertMatrix(proj_epsilon_matrix, inverse_proj_epsilon_matrix);
	result_of_mult= prod(inverse_proj_epsilon_matrix, diff_proj_i_j_vector);
	return inner_prod(diff_proj_i_j_vector, result_of_mult);
};

void SetNewWeights2NewTreeComparisons(
												std::vector<TreeComparison>& tree_comparison_vec,
												ublas::matrix<float, ublas::row_major> & epsilon_matrix,
												Parameters & pars,
												int index_of_begin) {
	int index_of_end = tree_comparison_vec.size();
	for(int i=index_of_begin; i < index_of_end; i++){

		float sum_denomination = 0;
		for(int j=0; j<index_of_begin; j++){
			sum_denomination += tree_comparison_vec[j].results.weight*
					get_gaussian_kernel(tree_comparison_vec,epsilon_matrix, i, j);
		};
		tree_comparison_vec[i].results.weight =
			1.0/(index_of_end*sum_denomination * 1.0);
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
																							pars.vars.gen);
	std::uniform_int_distribution<int> distribution(1,1000000);
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
	cout<<"cholesky_tria="<<L<<endl;  // end Cholesky_decompose

  int diff = pars.algo.parameter_population_num -
						 pars.algo.parameter_population_num_alpha;
	do{//shared(flag) #pragma omp parallel for
			int boundary = 2.0 * pars.algo.parameter_population_num;
  		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  		#pragma omp parallel for
			for(int i = 0; i < boundary; ++i) {
				SubInnerLoop(pars, tree_comparison_vec, epsilon_matrix, L);
			};
	} while( get_num_lines_in_comp_file(pars) <= diff);

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

void MainLoop(Parameters & pars,
							std::vector<TreeComparison>& tree_comparison_vec/*,
							mt19937 & gen*/)	{
	unsigned int iteratation_num = 1;
	pars.vars.p_acc = 1000;
	pars.vars.max_dist = std::max_element(tree_comparison_vec.begin(),
														 tree_comparison_vec.end(),
														 comparison_of_two_trees)->results.tree_dist;
	SetTreeComparisonVec2File(tree_comparison_vec, pars, 0);
	while ((pars.vars.max_dist > pars.algo.final_epsilon) &&
					(pars.vars.p_acc > pars.algo.p_acc_min))
	{
		std::sort(tree_comparison_vec.begin(), tree_comparison_vec.end(),
							comparison_of_two_trees);

	  WeightNormalization(tree_comparison_vec, 0,
												pars.algo.parameter_population_num_alpha);
		SetTreeComparisonVec2File(tree_comparison_vec, pars,
															iteratation_num);
		unsigned int m_size = pars.number_of_parameters;//matrix_size

		ublas::vector<float> weighted_mean_vector;
		weighted_mean_vector = GetWeightedMeanVec(tree_comparison_vec, m_size,
																			pars.algo.parameter_population_num_alpha);

		ublas::matrix<float, ublas::row_major> epsilon_matrix(m_size, m_size);
		epsilon_matrix = ublas::zero_matrix<float>(m_size, m_size);
		epsilon_matrix = get_twice_weighted_empirical_covariance(
																weighted_mean_vector,
																tree_comparison_vec,
																m_size,
																pars.algo.parameter_population_num_alpha);
		pars.vars.n_trials = 0;
		InnerLoop(pars, tree_comparison_vec, epsilon_matrix);

    float diff = tree_comparison_vec.size() -
                 pars.algo.parameter_population_num_alpha;
 		pars.vars.p_acc = (diff)/(pars.vars.n_trials*1.0);
		WeightNormalization(tree_comparison_vec,
												 pars.algo.parameter_population_num_alpha,
												 tree_comparison_vec.size());//w normalization
		WeightNormalization(tree_comparison_vec,0,
												 tree_comparison_vec.size());
		iteratation_num++;
		SetTreeComparisonVec2File(tree_comparison_vec, pars, iteratation_num);
	};
};




Parameters InitMainParameters(int argc, char *argv[]);

void DoesInitiationPart(Parameters & pars);
void GetPars(Parameters & pars);


int main (int argc, char *argv[]) {
	Parameters pars = InitMainParameters(argc, argv);
	DoesInitiationPart(pars);
	std::vector<TreeComparison> tree_comparison_vec;
  //get_tree_comparisons_from_file(m_par.NUM.c_str(),
  //array_of_tree_comparisons, m_par.parameter_population_num);
	GetTreeComparisonVecFromFile(tree_comparison_vec, pars);
  //init_weight(tree_comparison_vec, m_par.parameter_population_num);
	InitWeight(tree_comparison_vec, pars.algo.parameter_population_num);
#if defined(CHECK_REGIME)
	//set_tree_comparisons_to_file(array_of_tree_comparisons, pars, 100000);
	SetTreeComparisonVec2File(tree_comparison_vec, pars, 100000);
#endif // defined
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
};


