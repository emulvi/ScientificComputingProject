#ifndef HF_AUX_H
#define HF_AUX_H
#include <vector>
using namespace std;
//////Here are the functions required for HF

typedef std::vector<vector<double>> Real_Matrix;

double read_nuc_en(double nuc_en);

double read_T(int ao, Real_Matrix& T_int);

double read_S(int ao, Real_Matrix& S);

double read_v_int(int ao, Real_Matrix& v_int);

#endif
