
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
  float mut_rate;  // gama; // these are rates per daughter cell. Rates per diploid exome will be 2x higher (these values are given in the paper)
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
//        results.tree_dist = a.results.tree_dist;
//        results.delta_mut =  a.results.delta_mut;
//        results.num_mut =  a.results.num_mut;
//        results.weight = a.results.weight;
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
	float max_dist;
};

struct InfoSearch {
	float min_;
	float max_;
  float true_val;
};

struct SearchSpace{
  InfoSearch driver_adv;
  InfoSearch mut_rate;  // gama; // these are rates per daughter cell. Rates per diploid exome will be 2x higher (these values are given in the paper)
  InfoSearch driver_mut_rate; //driver_prob;
};

struct Parameters {
	string catalog_file_name;
	string par_file_name;
	string comparison_file_name;
	int seed; //Seed for initialization of Mersenne Twister pseudo-random number generator
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
void simulatate_and_write_to_file(const char* NUM, int i,  float max_dist,
																	int number_of_parameters, float driver_adv,
																	float mut_rate, float driver_mut_rate) {
    char command_line[256];
#if defined __linux
    sprintf(command_line,"./cancer_3d.exe %s 1 %d %f %f %f %f",NUM, i+10, driver_adv, driver_mut_rate, mut_rate , max_dist);
#else
    sprintf(command_line,"cancer_3d.exe %s 1 %d %f %f %f %f",NUM, i+10, driver_adv, driver_mut_rate, mut_rate, max_dist);
#endif
    cout <<"NUM="<<NUM<<endl;
    cout <<"command_line[256]="<<command_line<<endl;
    system (command_line);
    //cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!"<<j<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
};

int get_num_of_lines_in_file(const char* NUM)
{
    std::ifstream file_tree_comparison;
    char name_file[256];

    sprintf(name_file,"%s/file_tree_comparison.dat",NUM);
    file_tree_comparison.open(name_file);

    int number_of_lines = 0;
    std::string line;


    while (std::getline(file_tree_comparison, line))
        ++number_of_lines;
    return number_of_lines;
};

void get_tree_comparisons_from_file(const char* NUM,
                                    std::vector<TreeComparison> & array_of_tree_comparisons,
                                    int parameter_population_num) {
    std::ifstream file_tree_comparison;
    char name_file[256];

    sprintf(name_file,"%s/file_tree_comparison.dat",NUM);
    file_tree_comparison.open(name_file);

    if (file_tree_comparison.is_open()==false) {
			cout<<"err..."<<endl;
    } else {
			//do {
			for(int i = 0; i < parameter_population_num; i++) {
				TreeComparison A;

				file_tree_comparison >> A.par_point.driver_adv
														 >> A.par_point.driver_mut_rate
														 >> A.par_point.mut_rate
														 >> A.results.tree_dist
														 >> A.results.delta_mut
														 >> A.results.num_mut;
				if (!file_tree_comparison.eof()) array_of_tree_comparisons.push_back(A);
			};
			//} while (!file_tree_comparison.eof());
	};
};


void init_weight(std::vector <TreeComparison> & array_of_tree_comparisons,
								 int parameter_population_num) {
	for(int i=0; i < parameter_population_num; i++) {
		array_of_tree_comparisons[i].results.weight = 1.0/parameter_population_num;
	};
};

void set_tree_comparisons_to_file(const char* NUM,
																	std::vector<TreeComparison> &
																		array_of_tree_comparisons,
																	int number_of_iteration) {
	std::ofstream file_tree_comparison;
	char name_file[256];

	sprintf(name_file,"%s/file_tree_comparison_%d.dat", NUM, number_of_iteration);
	file_tree_comparison.open(name_file);

	if (file_tree_comparison.is_open()==false) cout<<"err..."<<endl;

	file_tree_comparison.precision(5);
	file_tree_comparison.setf(std::ios::fixed, std:: ios::floatfield);


	for(unsigned int i=0; i < array_of_tree_comparisons.size(); i++)
		file_tree_comparison<<"\t"<<
			array_of_tree_comparisons[i].par_point.driver_adv<<"\t"<<
			array_of_tree_comparisons[i].par_point.driver_mut_rate<<"\t"<<
			array_of_tree_comparisons[i].par_point.mut_rate<<"\t"<<
			array_of_tree_comparisons[i].results.tree_dist<<"\t"<<
			array_of_tree_comparisons[i].results.delta_mut<<"\t"<<
			array_of_tree_comparisons[i].results.num_mut<<"\t"<<
			array_of_tree_comparisons[i].results.weight<<endl;
};

void init_loop(Parameters pars) {
	char name_file[256];
	sprintf(name_file,"%s/file_tree_comparison.dat", pars.catalog_file_name.c_str());
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

	std::uniform_int_distribution<int> distribution(1,10000000);

	do {
			#pragma omp parallel for
			for(int i = 0; i <= pars.algo.parameter_population_num; ++i){
				simulatate_and_write_to_file(pars.catalog_file_name.c_str(),
																		 distribution(gen),
																		 1000,
																		 pars.number_of_parameters,
																		 driver_adv_gen(gen),
																		 mut_rate_gen(gen),
																		 driver_mut_rate_gen(gen));
			};
	} while (get_num_of_lines_in_file(pars.catalog_file_name.c_str())<=
						pars.algo.parameter_population_num);
};




void weight_normalization (std::vector<TreeComparison>& array_of_tree_comparisons,
													 int begin_index, int end_index) //parameter_population_num_alpha
{
    float weight_sum_alpha=0;
    for(int i = begin_index; i<end_index; i++)
			weight_sum_alpha+=array_of_tree_comparisons[i].results.weight;

   // cout<< "parameter_population_num_alpha="<<parameter_population_num_alpha;
   // cout<< "weight_sum_alpha="<<weight_sum_alpha;
    for(int i=begin_index; i<end_index; i++)
        array_of_tree_comparisons[i].results.weight=
        ((array_of_tree_comparisons[i].results.weight*1.0)/(1.0*weight_sum_alpha));
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

ublas::vector<float> get_weighted_mean_vector(std::vector<TreeComparison>& array_of_tree_comparisons,
                                              int size_of_matrix,
                                              int parameter_population_num_alpha) {
  ublas::vector<float> weighted_mean_vector(size_of_matrix);

  weighted_mean_vector = ublas::zero_vector<float> (size_of_matrix);
  weighted_mean_vector = std::accumulate(
                            array_of_tree_comparisons.begin(),
                            array_of_tree_comparisons.begin() +
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

  /////////////////////////

  for (int i=0; i<parameter_population_num_alpha; i++)	{
    ublas::vector<float> parameter_vector(size_matrix);

    parameter_vector(0) = array_of_tree_comparisons[i].par_point.driver_adv;
    parameter_vector(1) = array_of_tree_comparisons[i].par_point.mut_rate;
    parameter_vector(2) = array_of_tree_comparisons[i].par_point.driver_mut_rate;

    ublas::vector<float> diff_mean_vector(size_matrix);
    diff_mean_vector[0] = parameter_vector[0] - weighted_mean_vector[0];
    //change
    //if (abs(parameter_vector[1]-weighted_mean_vector[1])<0.0001)
    //diff_mean_vector[1]=0; else diff_mean_vector[1]=parameter_vector[1]-
	  //weighted_mean_vector[1];
    //if (abs(parameter_vector[2]-weighted_mean_vector[2])<0.0001)
    //diff_mean_vector[2]=0;
    //else diff_mean_vector[2]=parameter_vector[2]-
    //weighted_mean_vector[2];
    diff_mean_vector[1]=parameter_vector[1]-weighted_mean_vector[1];
    diff_mean_vector[2]=parameter_vector[2]-weighted_mean_vector[2];

    cout<<ublas::outer_prod(diff_mean_vector,diff_mean_vector)<<endl;
    epsilon_matrix += array_of_tree_comparisons[i].results.weight*
                       (ublas::outer_prod(diff_mean_vector,diff_mean_vector));
  };

  float sum_sq_wei = sum_of_square_weights(array_of_tree_comparisons,
																					 parameter_population_num_alpha);
  epsilon_matrix*=(2.0/(1-sum_sq_wei));


 // cout<<"\n"<<"weights... ";
 // for(int i=0;i<array_of_tree_comparisons.size();i++)
 // {j
 //       cout<<"a("<<i<<")="<<array_of_tree_comparisons.at(i).results.weight;
 // };
  //cout<<"\n"<<"size_vector= "<< weighted_mean_vector.size()<<endl;
 // std::vector<float> A=get_matrix_vector(epsilon_matrix, 0, size_matrix);
 // std::vector<std::vector<float>> AA=get_matrix(epsilon_matrix, size_matrix);

  cout<<"weighted_mean_vector="<<weighted_mean_vector<<endl;
  cout<<"epsilon_matrix="<<epsilon_matrix<<endl;


/*
  ublas::matrix<float, ublas::row_major> test_matrix(size_matrix, size_matrix);
  test_matrix=zero_m;


  cout<<"test_matrix=(zero)"<<endl;

  for (int i=0; i<size_matrix; i++) {
    for (int j = 0; j<size_matrix; j++) {
      cout<<test_matrix(i,j)<<"  ";
    }
    cout<<"\n";
  };

  parameter_vector(0) = array_of_tree_comparisons[0].parameters.driver_adv;
  parameter_vector(1) = array_of_tree_comparisons[0].parameters.gama;
  parameter_vector(2) = array_of_tree_comparisons[0].parameters.driver_prob;

    test_matrix+=(ublas::outer_prod(parameter_vector-weighted_mean_vector,parameter_vector-weighted_mean_vector));
  cout.precision(15);
  cout<<"\n"<<"parameter_vector="<<parameter_vector<<endl;
  cout<<"\n"<<"weighted_mean_vector="<<weighted_mean_vector<<endl;

  cout<<"parameter_vector-weighted_mean_vector="<<(parameter_vector-weighted_mean_vector)<<endl;


  cout<<"test_matrix=(multip)"<<endl;

  for (int i=0; i<size_matrix; i++) {
    for (int j = 0; j<size_matrix; j++) {
      cout<<test_matrix(i,j)<<"  ";
    }
    cout<<"\n";
  };

  test_matrix+=array_of_tree_comparisons[0].results.weight*(ublas::outer_prod(parameter_vector-weighted_mean_vector,parameter_vector-weighted_mean_vector));

  cout<<"test_matrix=(multip with weight)"<<endl;

  for (int i=0; i<size_matrix; i++) {
    for (int j = 0; j<size_matrix; j++) {
      cout<<test_matrix(i,j)<<"  ";
    }
    cout<<"\n";
  };
*/
 return epsilon_matrix;
};

int get_index_of_parameter_vector (
		std::vector<TreeComparison> array_of_tree_comparisons,
		int parameter_population_num_alpha, mt19937 & gen) {
	std::uniform_real_distribution<float> dis(0, 1);

	float uniform_coin=dis(gen);
	float sum = array_of_tree_comparisons[0].results.weight;
	int indicator = 0;

	//cout <<"sum="<<sum<<" indicator="<<indicator<<endl;
	while (sum < uniform_coin){
		sum += array_of_tree_comparisons[indicator].results.weight;
		indicator++;
	};

	//cout <<"sum_final="<<sum<<" indicator_final="<<indicator<<endl;
	return indicator;
};

ublas::vector<float>  get_multivariate_normal_distribution_sample(
												ublas::matrix<float, ublas::row_major> & L,
                        ublas::vector<float> & parameter_vector,
                        mt19937 & gen) {
	int size_matrix = L.size1();

	ublas::vector<float> new_parameter_vector(size_matrix);

	//search acceptable parameters
	//do {
	std::normal_distribution<float> n_d(0,1);

	new_parameter_vector[0] = n_d(gen);
	new_parameter_vector[1] = n_d(gen);
	new_parameter_vector[2] = n_d(gen);

	new_parameter_vector = prod( new_parameter_vector, L);
	new_parameter_vector += parameter_vector;

	if  (new_parameter_vector[0]<0) new_parameter_vector[0]= -1*new_parameter_vector[0];
	if  (new_parameter_vector[1]<0) new_parameter_vector[1]= -1*new_parameter_vector[1];
	if  (new_parameter_vector[2]<0) new_parameter_vector[2]= -1*new_parameter_vector[2];
		//change
		//if ( (new_parameter_vector[0]>=0) && (new_parameter_vector[1]>=0) && (new_parameter_vector[2]>=0) &&
		//     (new_parameter_vector[0]<=1) && (new_parameter_vector[1]<=1) && (new_parameter_vector[2]<=1)
		//    ) new_parmeters_in_domain = true;

	//} while (new_parmeters_in_domain == false);
	//--------------------------------------------------
	//change
	//if (new_parmeters_in_domain == true) return new_parameter_vector;

	//return parameter_vector;
	return new_parameter_vector;
};

float get_gaussian_kernel(
				std::vector<TreeComparison>& array_of_tree_comparisons,
				ublas::matrix<float, ublas::row_major> & epsilon_matrix,
				int i, int j) {
	int number_of_parameters = 1;
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

	 //cout<<"inverse_proj_epsilon_matrix=" <<inverse_proj_epsilon_matrix<<endl;

	return inner_prod(diff_proj_i_j_vector, result_of_mult);
};

void set_new_weights_to_new_tree_comparisons(int index_of_begin,
				std::vector<TreeComparison>& array_of_tree_comparisons,
				ublas::matrix<float, ublas::row_major> & epsilon_matrix,
				mt19937 & gen ) {
	int index_of_end = array_of_tree_comparisons.size();
	for(int i=index_of_begin; i < index_of_end; i++){

		float sum_denomination = 0;
		for(int j=0; j<index_of_begin; j++){
			sum_denomination += array_of_tree_comparisons[j].results.weight*
					get_gaussian_kernel(array_of_tree_comparisons,epsilon_matrix, i, j);
		};
		array_of_tree_comparisons[i].results.weight =
			1.0/(index_of_end*sum_denomination*1.0);

		//cout<<"for i="<<i<<" the sum_denomination="<<sum_denomination
		//<<"array_of_tree_comparisons.results.weight "<<array_of_tree_comparisons[i].results.weight<<endl;

		//cout<< array_of_tree_comparisons<<endl;
	};
};



void sub_inner_loop(int parameter_population_num_alpha,
										std::vector<TreeComparison> & array_of_tree_comparisons,
                    ublas::matrix<float, ublas::row_major> & epsilon_matrix,
                    const char* NUM, mt19937 & gen , int & N_trials,
                    ublas::matrix<float, ublas::row_major> Cholesky_matrix,
                    float max_dist, int  number_of_parameters) {
	int matrix_size = epsilon_matrix.size1();

	int index = get_index_of_parameter_vector (array_of_tree_comparisons,
																						 parameter_population_num_alpha,
																						 gen);

	ublas::vector<float> parameter_vector(matrix_size);
	//cout<<"index of randomly taken parameters with good tree dist ="<<index<<endl;

	parameter_vector[0] = array_of_tree_comparisons[index].par_point.driver_adv;
	parameter_vector[1] = array_of_tree_comparisons[index].par_point.mut_rate;
	parameter_vector[2] = array_of_tree_comparisons[index].par_point.driver_mut_rate;

	//cout <<"parameter_vector="<<parameter_vector<<endl;

	parameter_vector = get_multivariate_normal_distribution_sample(Cholesky_matrix,
																																parameter_vector,
																																gen);

	std::uniform_int_distribution<int> distribution(1,1000000);
	TreeComparison tree_com;
	tree_com.par_point.driver_adv = parameter_vector [0];
	tree_com.par_point.mut_rate = parameter_vector [1];
	tree_com.par_point.driver_mut_rate = parameter_vector [2];

	simulatate_and_write_to_file(NUM, distribution(gen), max_dist,
															 number_of_parameters,
															 tree_com.par_point.driver_adv,
															 tree_com.par_point.mut_rate,
															 tree_com.par_point.driver_mut_rate);

	N_trials++;
};

void inner_loop(int parameter_population_num_alpha,
								std::vector<TreeComparison> & array_of_tree_comparisons,
                ublas::matrix<float, ublas::row_major> & epsilon_matrix,
                const char* NUM, mt19937 & gen,int & N_trials,
                int  number_of_parameters) {
	int matrix_size = epsilon_matrix.size1();
	float max_dist = array_of_tree_comparisons[
												 array_of_tree_comparisons.size() -
	                       parameter_population_num_alpha-1
																						].results.tree_dist;

 // cout<<"\n"<<"array_of_tree_comparisons.size()-parameter_population_num_alpha-1="<<array_of_tree_comparisons.size()-parameter_population_num_alpha-1
 // <<"  max_dist="<<max_dist <<endl;

	//////////////////////////////////////////////////////init file for new loop//////////////////////////////////////////////////////
	char name_file[256];
	sprintf(name_file,"%s/file_tree_comparison.dat",NUM);
	remove(name_file);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//volatile bool flag = false;



	///////////////////////////////////////////////////////// cholesky_decompose ////////////////////////////////////////////////////////////////
	ublas::matrix<float, ublas::row_major> L (matrix_size, matrix_size);
	L=ublas::zero_matrix<float>(matrix_size, matrix_size);
	cholesky_decompose(epsilon_matrix, L);
	cout<<"cholesky_tria="<<L<<endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int parameter_population_num = array_of_tree_comparisons.size();



	do{//shared(flag) #pragma omp parallel for
			int boundary = 2.0 * parameter_population_num;
			#pragma omp parallel for
			for(int i = 0; i < boundary; ++i) //
			{
			//    if(flag) continue;
					sub_inner_loop(parameter_population_num_alpha, array_of_tree_comparisons,
					epsilon_matrix, NUM, gen , N_trials, L, max_dist,  number_of_parameters);
			};

	} while (get_num_of_lines_in_file(NUM)<=
						(parameter_population_num - parameter_population_num_alpha));

	std::vector<TreeComparison> additional_array_of_tree_comparisons;
	get_tree_comparisons_from_file(NUM, additional_array_of_tree_comparisons,
										(parameter_population_num - parameter_population_num_alpha));

	//send results to main array array_of_tree_comparisons from additional array
	for(int i=parameter_population_num_alpha; i<parameter_population_num;i++)
		array_of_tree_comparisons[i]=additional_array_of_tree_comparisons
																	[i-parameter_population_num_alpha];


	set_new_weights_to_new_tree_comparisons(parameter_population_num_alpha,
																					array_of_tree_comparisons,
																					epsilon_matrix, gen);
}

void main_loop( Parameters pars,
								std::vector<TreeComparison>& array_of_tree_comparisons,
                mt19937 & gen)	{
	int number_of_iteration = 1;
	float max_dist;
	float p_acc=1000;

	max_dist=std::max_element(array_of_tree_comparisons.begin(),
														 array_of_tree_comparisons.end(),
														 comparison_of_two_trees)->results.tree_dist;
	set_tree_comparisons_to_file(pars.catalog_file_name.c_str(),
																array_of_tree_comparisons,0);
	while ((max_dist > pars.algo.final_epsilon) && (p_acc > pars.algo.p_acc_min)){
		std::sort(array_of_tree_comparisons.begin(),
								array_of_tree_comparisons.end(),
								comparison_of_two_trees);

		weight_normalization(array_of_tree_comparisons, 0,
												pars.algo.parameter_population_num_alpha);//w normalization
		set_tree_comparisons_to_file( pars.catalog_file_name.c_str(),
														 	    array_of_tree_comparisons,number_of_iteration);

		//////////////////////////////////
		unsigned int matrix_size=3;
		ublas::matrix<float, ublas::row_major> epsilon_matrix(matrix_size, matrix_size);
		epsilon_matrix = ublas::zero_matrix<float>(matrix_size, matrix_size);
		ublas::vector<float> weighted_mean_vector;

		weighted_mean_vector = get_weighted_mean_vector(array_of_tree_comparisons,
																			matrix_size,
																			pars.algo.parameter_population_num_alpha);

		epsilon_matrix = get_twice_weighted_empirical_covariance(
																weighted_mean_vector,
																array_of_tree_comparisons,
																matrix_size,
																pars.algo.parameter_population_num_alpha);

		int N_trials = 0;
		inner_loop(pars.algo.parameter_population_num_alpha,
							 array_of_tree_comparisons, epsilon_matrix,
							 pars.catalog_file_name.c_str(), gen, N_trials,
							 pars.number_of_parameters);
    float diff = array_of_tree_comparisons.size() -
                 pars.algo.parameter_population_num_alpha;
 		p_acc = (diff)/(N_trials*1.0);
		weight_normalization(array_of_tree_comparisons,
												 pars.algo.parameter_population_num_alpha,
												 array_of_tree_comparisons.size());//w normalization
		weight_normalization(array_of_tree_comparisons,0,
												 array_of_tree_comparisons.size());

				//cout<<"par vec="<<parameter_vector <<endl;
				//boost::multivariate_normal_distribution<float> m_n_d(L,weighted_mean_vector);
				//boost::m

				//multivariate_normal_distribution  (const matrix_type& cholesky, const vector_type& mean)
				//cout<<endl<<"cholesky_decompose of epsilon_matrix="<<L<<endl;

		//    } do;

		//};






		number_of_iteration++;
		set_tree_comparisons_to_file(pars.catalog_file_name.c_str(),
																 array_of_tree_comparisons,
																 number_of_iteration);
		//cout<<"\n"<<" p_acc="<<p_acc<<"; max_dist="<<std::max_element(array_of_tree_comparisons.begin(), array_of_tree_comparisons.end(), comparison_of_two_trees)->results.tree_dist<<endl;
   };
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
  const int expected_arg_num = 2; // the number of arguments in command line
  if (argc == expected_arg_num) {
		pars.par_file_name = argv[1];
		pars.seed = atof(argv[2]);
		return pars;
	} else {
		err("MISTAKE: there are 2 arguments: first - file name; second - seed ");
	};
	return pars;
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
		("algorithm_parameters.parameter_population_num");
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
		("search_space.mut_rate.min");
  pars.search_space.mut_rate.max_ =  pt.get<float>
		("search_space.mut_rate.max");
  pars.search_space.mut_rate.true_val = pt.get<float>
		("search_space.mut_rate.true");

  pars.search_space.driver_mut_rate.min_ = pt.get<float>
		("search_space.driver_mut_rate.min");
  pars.search_space.driver_mut_rate.max_ = pt.get<float>
		("search_space.driver_mut_rate.max");
  pars.search_space.driver_mut_rate.true_val = pt.get<float>
    ("search_space.driver_mut_rate.true");
};

void DoesInitiationPart(Parameters & pars) {
	string full_par_file_name = pars.par_file_name;

	if (FileExists(full_par_file_name.c_str())) {
		GetPars(pars);
	} else {
		cout << "file with the name "<<  pars.par_file_name <<
						" doesn't exist"<< endl;
	};
};


int main (int argc, char *argv[]) {
	//Parameters pars = InitMainParameters(argc, argv);

//ofstream matrix_file;
// string full_comparison_file_name = catalog_name + "/" +
//                                   file_name+"Mutations" +".info";
//matrix_file.open(full_comparison_file_name.c_str());

////////#if defined __linux
////////    if ( boost::filesystem::exists("./model_info.info" ) ){
////////#else
////////    if ( boost::filesystem::exists("model_info.info" ) ){
////////#endif
//#if defined __linux
//   std::ifstream infile("./model_info.info");
//#else
//   std::ifstream infile("model_info.info");
//#endif

//if (infile.good()){
//    //initiation
//
//
//    init_prameters m_par = init_par();
//    mt19937 gen; //Initialize Mersenne Twister pseudo-random number generator
//    gen.seed(m_par.seed);
//    init_loop(m_par); //simulate parameters and results of parameter_population_num and put it in file_tree_comparison.dat
//
//
//    std::vector<tree_comparison> array_of_tree_comparisons;
//    get_tree_comparisons_from_file(m_par.NUM.c_str(), array_of_tree_comparisons, m_par.parameter_population_num);
//
//    init_weight(array_of_tree_comparisons, m_par.parameter_population_num);
//
//
//#if defined(CHECK_REGIME)
//    set_tree_comparisons_to_file(m_par.NUM.c_str(), array_of_tree_comparisons,100000);
//#endif // defined
//
//
//    main_loop(m_par, array_of_tree_comparisons, gen);
//    //int number_of_lines = get_num_of_lines_in_file(NUM);
//    //std::cout << "Number of lines in text file: " << number_of_lines;
//
//} else {
//    std::cout << "Can't find file model_info.info! " << std::endl;
//};
  cout<<10<<endl;
  return 0;
};





