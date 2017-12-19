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


//====
// 2D
//====

   struct M4l_vs_Dkinbkg
   {
      TString x_label = "m_{4#font[12]{l}} (GeV)";
      TString y_label = "D_{bkg}^{kin}";
      TString cut_label = "";
      int x_bins= 100;
      float x_min = 50;
      float x_max = 850;
      int y_bins= 20;
      float y_min = 0;
      float y_max = 1;
   };

   struct MZ2_vs_Dkinbkg
   {
      TString x_label = "m_{Z_{2}} (GeV)";
      TString y_label = "D_{bkg}^{kin}";
      TString cut_label = "";
      int x_bins= 75;
      float x_min = 0;
      float x_max = 150;
      int y_bins= 20;
      float y_min = 0;
      float y_max = 1;
   };
   
      struct D2jVbfHjj_vs_D2jqg
   {
      TString x_label = "D_{VBF-2j}^{ME}";
      TString y_label = "D_{2jets}^{q/g}";
      TString cut_label = "";
      int x_bins= 51;
      float x_min = 0;
      float x_max = 1.02;
      int y_bins= 25;
      float y_min = -0.04;
      float y_max = 0.96;
   };

};
#endif
