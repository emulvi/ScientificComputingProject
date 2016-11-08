#ifndef HF_AUX_H
#define HF_AUX_H
#include <iostream>
#include <iomanip> 
#include <cstdlib> 
#include <vector>
#include <cmath>

using namespace std;
//////Here are the functions required for HF

typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

void read_nuc_en(double nuc_en);

void read_T(int ao, Real_Matrix& T_int);

void read_S(int ao, Real_Matrix& S);

void read_v_nuc(int ao, Real_Matrix& v_nuc);

void build_H_core(int ao, Real_Matrix& v_nuc, Real_Matrix& T_int, Real_Matrix& H_core);

void read_v_int(int ao, Real_4dMatrix& v_int);

#endif
