#ifndef DEBUG_FUNCTIONS_H
#define DEBUG_FUNCTIONS_H

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>


namespace ublas = boost::numeric::ublas;

float get_element(ublas::matrix<float, ublas::row_major> u_matrix, unsigned int i, unsigned int j);
std::vector<float> get_matrix_vector(ublas::matrix<float, ublas::row_major> u_matrix, unsigned int row, unsigned int matrix_size);
std::vector<std::vector<float>> get_matrix(ublas::matrix<float, ublas::row_major> u_matrix,  unsigned int matrix_size);

std::vector<float> get_vector(const ublas::vector<float> u_vec);




#endif // DEBUG_FUNCTIONS_H
