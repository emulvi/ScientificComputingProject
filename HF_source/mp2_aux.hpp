//////Here are the functions required for mp2
//
#ifndef MP2_AUX_H
#define MP2_AUX_H
#include <iostream>
#include <iomanip> 
#include <cstdlib> 
#include <vector>
#include <cmath>

using namespace std;
//////Here are the functions required for HF

typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;



void transform_v_int(int ao, Matrix& C, Matrix& v_int, Matrix& v_int_mo, Matrix& Xmat, Matrix& evecs);

void transform_v_int_2(int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo);

void transform_v_int_CD(int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo);

double calculate_E_mp2(int ao, int occ, Matrix& evals, Matrix& v_int_mo);

#endif
