#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include "hf_aux.h"

//////Here are the functions required for HF

double read_nuc_en(double nuc_en)
{

   std::ifstream nuc_ener("enuc.dat");
   while(!nuc_ener.eof()){
      nuc_ener >> nuc_en;
   }

//   std::cout << "Nuclear energy is: " << nuc_en << std::endl;

   return nuc_en;

}


double read_T(int ao){

   double T_int[ao-1][ao-1];
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
      return T_int[ao-1][ao-1];
   }


}

double read_S(int ao){

   double S[ao-1][ao-1];
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
   return S[ao-1][ao-1];


}

double read_v_int(int ao){

   double v_int[ao-1][ao-1];
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
   return v_int[ao-1][ao-1];

}
