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



void transform_v_int(int ao, Matrix& C, Real_4dMatrix& v_int, Real_4dMatrix& v_int_mo);

double calculate_E_mp2(int ao, int occ, Matrix& evals, Real_4dMatrix& v_int_mo);

#endif
