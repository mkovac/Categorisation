// Include classes
#include "Variables.h"

using namespace std;

// Constructor
//=======================
Variables::Variables() {

   vec_var_pairs.push_back(variable_pairs());
   vec_var_pairs[0].x_label = "m_{4#font[12]{l}} (GeV)";
   vec_var_pairs[0].y_label = "D_{bkg}^{kin}";
   vec_var_pairs[0].x_bins = 100;
   vec_var_pairs[0].x_min  = 50;
   vec_var_pairs[0].x_max  = 850;
   vec_var_pairs[0].y_bins = 20;
   vec_var_pairs[0].y_min  = 0;
   vec_var_pairs[0].y_max  = 1;
   
   vec_var_pairs.push_back(variable_pairs());
   vec_var_pairs[1].x_label = "m_{Z_{2}} (GeV)";
   vec_var_pairs[1].y_label = "D_{bkg}^{kin}";
   vec_var_pairs[1].x_bins = 75;
   vec_var_pairs[1].x_min  = 0;
   vec_var_pairs[1].x_max  = 150;
   vec_var_pairs[1].y_bins = 20;
   vec_var_pairs[1].y_min  = 0;
   vec_var_pairs[1].y_max  = 1;
   
   vec_var_pairs.push_back(variable_pairs());
   vec_var_pairs[2].x_label = "D_{VBF-2j}^{ME}";
   vec_var_pairs[2].y_label = "D_{2jets}^{q/g}";
   vec_var_pairs[2].x_bins = 51;
   vec_var_pairs[2].x_min  = 0;
   vec_var_pairs[2].x_max  = 1.02;
   vec_var_pairs[2].y_bins = 25;
   vec_var_pairs[2].y_min  = -0.04;
   vec_var_pairs[2].y_max  = 0.96;

}
//=======================


//========================
Variables::~Variables() {}
//========================
