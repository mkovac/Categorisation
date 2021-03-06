#ifndef ROC_h
#define ROC_h

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TString.h"


using namespace std;

class ROC
{
   
public:
   ROC();
   ~ROC();
   void Prepare( int, float , int , int , float , bool );
   void Clean();


   struct ROCs
   {
      int ROC_num;
      float var_value;
      int sig_proc;
      int bkg_proc;
      float jet_cut;
      bool cut;
   };
   
   vector<ROCs> vec_ROCs;
   ROCs temp_ROC;
};
#endif
