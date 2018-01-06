#ifndef Variables_h
#define Variables_h

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TString.h"


using namespace std;

class Variables
{
   
public:
   Variables();
   ~Variables();
   void PrepareVarPair( TString, TString, TString, int, int, int, int, int, int );
   void PrepareVar( TString, TString, int, int, int, int);
   
   struct pt
   {
      TString var_X_label = "p_{T}^{#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 1";
      TString var_cut_label = "";
      int var_N_bin   = 80;
      float var_min   = 0;
      float var_max   = 80;
      bool var_log_x  = 0;
      bool var_log_y  = 0;
      int var_CMS_pos = 33;
      int varLegPos   = 33;
   };
   
   struct eta
   {
      TString var_X_label = "|#eta|^{#font[12]{l}}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      int var_N_bin = 100;
      float var_min = 0;
      float var_max = 5;
      bool var_log_x = 0;
      bool var_log_y = 0;
      int var_CMS_pos = 11;
      int varLegPos = 33;
   };
   
   struct n_ele
   {
      TString var_X_label = "# of reco electrons";
      TString var_Y_label = "Events / 1";
      TString var_cut_label = "";
      int var_N_bin = 10;
      float var_min = 0;
      float var_max = 10;
      bool var_log_x = 0;
      bool var_log_y = 0;
      int var_CMS_pos = 11;
      int varLegPos = 33;
   };
   
   struct n_mu
   {
      TString var_X_label = "# of reco muons";
      TString var_Y_label = "Events / 1";
      TString var_cut_label = "";
      int var_N_bin = 10;
      float var_min = 0;
      float var_max = 10;
      bool var_log_x = 0;
      bool var_log_y = 0;
      int var_CMS_pos = 11;
      int varLegPos = 33;
   };
   
   struct n_lep
   {
      TString var_X_label = "# of reco leptons";
      TString var_Y_label = "Events / 1";
      TString var_cut_label = "";
      int var_N_bin = 10;
      float var_min = 0;
      float var_max = 10;
      bool var_log_x = 0;
      bool var_log_y = 0;
      int var_CMS_pos = 11;
      int varLegPos = 33;
   };


//===============
// Variable pairs
//===============
  
   struct variable_pairs
   {
      TString name;
      TString x_label;
      TString y_label;
      TString cut_label;
      int x_bins;
      float x_min;
      float x_max;
      int y_bins;
      float y_min;
      float y_max;
   };
   
   vector<variable_pairs> vec_var_pair;
   variable_pairs temp_var_pair;


//==========
// Variables
//==========
  
   struct variables
   {
      TString name;
      TString label;
      TString cut_label;
      int n_bins;
      float x_min;
      float x_max;
      float x_max_roc;
   };

   vector<variables> vec_var;
   variables temp_variable;
};
#endif
