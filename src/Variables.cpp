// Include classes
#include "Variables.h"

using namespace std;

// Constructor
//======================
Variables::Variables(){}
//======================


//=======================
Variables::~Variables(){}
//=======================


//==================================================================================================================================================
void Variables::PrepareVarPair( TString name, TString x_label, TString y_label, int x_bins, int x_min, int x_max, int y_bins, int y_min, int y_max )
{
   temp_var_pair.name    = name;
   temp_var_pair.x_label = x_label;
   temp_var_pair.y_label = y_label;
   temp_var_pair.x_bins  = x_bins;
   temp_var_pair.x_min   = x_min;
   temp_var_pair.x_max   = x_max;
   temp_var_pair.y_bins  = y_bins;
   temp_var_pair.y_min   = y_min;
   temp_var_pair.y_max   = y_max;
   
   vec_var_pair.push_back(temp_var_pair);
}
//==================================================================================================================================================



//=========================================================================================
void Variables::PrepareVar( TString name, TString label, int n_bins, int x_min, int x_max )
{
   temp_variable.name   = name;
   temp_variable.label  = label;
   temp_variable.n_bins = n_bins;
   temp_variable.x_min  = x_min;
   temp_variable.x_max  = x_max;
   
   vec_var.push_back(temp_variable);
}
//=========================================================================================

