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
   Variables ();
   ~Variables();

   
//====
// pt
//====
   
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
   
};
#endif
