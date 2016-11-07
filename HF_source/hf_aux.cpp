#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include "hf_aux.h"
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;
//////Here are the functions required for HF

typedef std::vector<vector<double>> Real_Matrix;

void read_nuc_en(double nuc_en)
{
   std::ifstream nuc_ener("enuc.dat");
   while(!nuc_ener.eof()){
      nuc_ener >> nuc_en;
   }

//   std::cout << "Nuclear energy is: " << nuc_en << std::endl;

   return;

}


void read_T(int ao, Real_Matrix& T_int){
   
  //Real_Matrix T_int(ao-1,vector<double>(ao-1,0.0));
   //double T_int[ao-1][ao-1];
   double val;
   int i;
   int j;

   std::ifstream kin_en;
   kin_en.open("T.dat");
   for(int k=0;k<(ao*(ao+1))/2;k++){
      kin_en >> i;
      kin_en >> j;
      kin_en >> val;
      T_int[i-1][j-1] = val;
      std::cout << i << " " << j << " " << T_int[i-1][j-1] << std::endl;
      std::cout << typeid(T_int[ao-1][ao-1]) << std::endl;
   }
   return;
}

void read_S(int ao, Real_Matrix& S){

   Real_Matrix S(ao-1,vector<double>(ao-1,0.0));
   //double S[ao-1][ao-1];
   double val;
   int i;
   int j;

   std::ifstream overlap;
   overlap.open("overlap.dat");
   for(int k=0;k<(ao*(ao+1))/2;k++){
      overlap >> i;
      overlap >> j;
      overlap >> val;
      S[i-1][j-1] = val;
      std::cout << i << " " << j << " " << S[i-1][j-1] << std::endl;
   }
   return;
}

void read_v_int(int ao, Real_Matrix& v_int){

   Real_Matrix v_int(ao-1,vector<double>(ao-1,0.0));
   //double v_int[ao-1][ao-1];
   double val;
   int i;
   int j;


   std::ifstream interaction;
   interaction.open("v_int.dat");
   for(int k=0;k<(ao*(ao+1))/2;k++){
      interaction >> i;
      interaction >> j;
      interaction >> val;
      v_int[i-1][j-1] = val;
      std::cout << i << " " << j << " " << v_int[i-1][j-1] << std::endl;
   }
   return;

}
