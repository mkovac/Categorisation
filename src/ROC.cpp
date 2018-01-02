// Include classes
#include "ROC.h"

using namespace std;

// Constructor
//==========
ROC::ROC(){}
//==========


//===========
ROC::~ROC(){}
//===========


//==================================================================================================
void ROC::Prepare( TString name, float var_value, int sig_proc, int bkg_proc, float jet_cut, bool cut )
{

   temp_ROC.name      = name;
   temp_ROC.var_value = var_value;
   temp_ROC.sig_proc  = sig_proc;
   temp_ROC.bkg_proc  = bkg_proc;
   temp_ROC.jet_cut   = jet_cut;
   temp_ROC.cut       = cut;

   vec_ROCs.push_back(temp_ROC);
}
//==================================================================================================



//===============
void ROC::Clean()
{
   vec_ROCs.clear();
}
//===============
