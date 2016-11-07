#ifndef HF_AUX_H
#define HF_AUX_H
#include <iostream>
#include <iomanip> 
#include <cstdlib> 
#include <vector>

using namespace std;
//////Here are the functions required for HF

typedef std::vector<vector<double>> Real_Matrix;

void read_nuc_en(double nuc_en);

void read_T(int ao, Real_Matrix T_int);

void read_S(int ao, Real_Matrix S);

void read_v_int(int ao, Real_Matrix v_int);

#endif
