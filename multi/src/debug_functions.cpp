#include "debug_functions.h"

float get_element(ublas::matrix<float, ublas::row_major> u_matrix, unsigned int i, unsigned int j)
{
    return u_matrix(i,j);
}

std::vector<float> get_vector(const ublas::vector<float> u_vec)
{
        std::vector<float> vec;

        for(unsigned int i=0; i<u_vec.size(); i++)
        {
          vec.push_back(u_vec(i));
        };
    return vec;

}

std::vector<float> get_matrix_vector(ublas::matrix<float, ublas::row_major> u_matrix, unsigned int row, unsigned int matrix_size)
{
        std::vector<float> matrix_vector;
        for(unsigned int i=0; i<matrix_size; i++)
        {
          matrix_vector.push_back(u_matrix(row,i));
        };
    return matrix_vector;
}

std::vector<std::vector<float>> get_matrix(ublas::matrix<float, ublas::row_major> u_matrix,  unsigned int matrix_size)
{
        std::vector<std::vector<float>> this_matrix;
        for(unsigned int i=0; i<matrix_size; i++)
        {
          std::vector<float> matrix_vector;
          this_matrix.push_back(matrix_vector);
          for(unsigned int j=0; j<matrix_size; j++) this_matrix.at(i).push_back(u_matrix(i,j));
        };
    return this_matrix;
}
