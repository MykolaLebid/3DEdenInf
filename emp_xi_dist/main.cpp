#include <iostream>
//parallel calculations
#include <omp.h>

#include <random>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element
#include <time.h>       /* time */

using namespace std;

int main()
{
		const int number_of_trials = 2000;
		int max_dimention = 100000;

		std::random_device rd_1;
    std::mt19937 mt1(rd_1()), mt2(time(NULL));
		std::uniform_int_distribution<int> rnd(1,1000000000);

		vector <int> seed_vec_1(number_of_trials);
		vector <int> seed_vec_2(number_of_trials);

		for (auto & x : seed_vec_1) x = rnd(mt1);
		for (auto & x : seed_vec_2) x = rnd(mt2);

//		for (auto & x : seed_vec_1) std::cout<<x<<std::endl;
//		std::cout << std::endl;
//		for (auto & x : seed_vec_2) std::cout<<x<<std::endl;

		for(int dim = 1; dim <= max_dimention; dim++) {
				vector<int> sample_numbers_1(dim, 0);
				vector<int> sample_numbers_2(dim, 0);

				vector<float> xi_sq_vec(number_of_trials,0);
				//#pragma omp parallel for
				for(int j = 0; j < number_of_trials; j++) {
						vector<float> weights_d_1(dim, 0);
						vector<float> weights_d_2(dim, 0);

						std::uniform_int_distribution <int> distribution_d_1(0,10),
																								distribution_d_2(0,10);
						// Setup the random vectors
						float sum_of_vec_1_elemets = 0,
								  sum_of_vec_2_elemets = 0; // scaling factors
						std::mt19937 mt_1(seed_vec_1.at(j)), mt_2(seed_vec_2.at(j));

						for (int i = 0; i < dim; i++) {
							sample_numbers_1.at(i) =  distribution_d_1(mt_1);
							sample_numbers_2.at(i) =  distribution_d_2(mt_2);

//							std::cout<<sample_numbers_1.at(i)<<"; "
//											 <<sample_numbers_2.at(i)<<std::endl;

							sum_of_vec_1_elemets += sample_numbers_1.at(i);
							sum_of_vec_2_elemets += sample_numbers_2.at(i);
						}

//						std::cout<< "sum_of_vec_1_elemts = " << sum_of_vec_1_elemets<<
//											  ";sum_of_vec_2_elemts = " << sum_of_vec_2_elemets<<std::endl;
				    // scale mass to 1

						for (int i = 0; i < dim; i++) {
							weights_d_1.at(i) = sample_numbers_1.at(i)*1.0/ sum_of_vec_1_elemets;

							weights_d_2.at(i) = sample_numbers_2.at(i)*1.0/ sum_of_vec_2_elemets;
							//std::cout<<weights_d_2.at(i);
						}
						//std::cout <<std::endl;
						float init = 0;

						float xi_sq = 0;
						for (int i = 0; i < dim; i++) {
							 float v_1 = weights_d_1.at(i);
							 float v_2 = weights_d_2.at(i);
							 if ((v_1 != 0) || (v_2 != 0)) {
									xi_sq += pow(v_1 - v_2, 2.0)/(v_1/2 + v_2/2);
							 };
						};
						xi_sq_vec.at(j) = xi_sq;
//						std::cout << xi_sq << endl;

				};
				auto max_e = std::max_element(xi_sq_vec.begin(), xi_sq_vec.end());
				auto min_e = std::min_element(xi_sq_vec.begin(), xi_sq_vec.end());
				std::cout<<"dim="<<dim<<"; ["<<(*min_e)<<";"<<(*max_e)<<"]" << std::endl;
				//TODO find max and min of xi_sq_vec and print it

		};
    return 0;
}
