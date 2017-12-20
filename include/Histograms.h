#ifndef Histograms_h
#define Histograms_h

// C++
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <tuple>

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TIterator.h"
#include "TROOT.h"

// Include classes
#include "CMS_lumi.h"
#include "Variables.h"
#include "Cosmetics.h"
#include "Counters.h"


using namespace std;

//const int num_of_production_modes    = Settings::num_of_production_modes;
//const int num_of_processes           = Settings::num_of_processes;
//const int num_of_processes_yields    = Settings::num_of_processes_yields;
//const int num_of_final_states        = Settings::num_of_final_states;
//const int num_of_categories          = Settings::num_of_categories;
//const int num_of_1D_plot_names       = Settings::num_of_1D_plot_names;
//const int num_of_2D_plot_names       = Settings::num_of_2D_plot_names;
//const int num_of_2D_error_plot_names = Settings::num_of_2D_error_plot_names;


class Histograms
{
   
public:
   Histograms(float);
   ~Histograms();
   
   void FillPt( float, float, int, int, int, int, int );
   void FillPtReco( float, float, int, int, int, int, int );
   void FillAbsEta( float, float, int, int, int, int, int );
   void FillAbsEtaReco( float, float, int, int, int, int, int );
   void FillPtEtaH( float, float, float, int );
   
   void FillRecoEleN( int, float, int, int, int, int );
   void FillRecoMuN( int, float, int, int, int, int );
   void FillRecoLepN( int, float, int, int, int, int );
   
   void FillVariables( float, float, int, int, int, int );
   void FillVariablePairs( float, float, float, int, int, int, int );
   void FillVariablePairsDecay( float, float, float, int, int, int, int, int );

   void SaveHistograms( TString );
   

private:
      
   float lumi_;
   TString histo_name_, histo_label_;
   vector<TString> s_process_, s_sort_, s_gen_ch_, s_reco_ch_, s_variable_, s_variable_pair_, s_variable_decay_;
   
   vector<tuple<int, float, float>> var_bins_min_max_;

   tuple<int, float, float> tpl_;
   
   Variables var_pairs;
   
//==============
// 1D histograms
//==============
   TH1F *h_pt_gen_H_lep_in_eta_acc_[Counters::num_of_processes][Counters::num_of_gen_channels][Counters::num_of_sorted_objects];
   TH1F *h_pt_reco_bc_in_sig_reg_and_pass_triger_[Counters::num_of_processes][Counters::num_of_gen_channels][Counters::num_of_sorted_objects];
   TH1F *h_eta_gen_H_lep_in_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels][Counters::num_of_sorted_objects];
   TH1F *h_eta_reco_bc_in_sig_reg_and_pass_triger_[Counters::num_of_processes][Counters::num_of_gen_channels][Counters::num_of_sorted_objects];

   TH1F *h_num_reco_H_ele_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   TH1F *h_num_reco_H_mu_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   TH1F *h_num_reco_H_lep_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   
   TH1F *h_gen_H_pt_[Counters::num_of_processes];
   TH1F *h_gen_H_eta_[Counters::num_of_processes];
   TH2F *h_gen_H_eta_vs_pt_[Counters::num_of_processes];
   
   TH1F* h_bc_in_sig_reg_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH1F* h_bc_in_sig_reg_match_H_leps_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH1F* h_bc_in_sig_reg_match_all_leps_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH1F* h_bc_in_sig_reg_match_WH_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH1F* h_bc_in_sig_reg_match_ZH_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH1F* h_bc_in_sig_reg_match_ttH_[Counters::num_of_variables][Counters::num_of_processes][Counters::num_of_reco_channels];

   TH2F* h_bc_in_sig_reg_2D_[Counters::num_of_variable_pairs][Counters::num_of_processes][Counters::num_of_reco_channels];
   TH2F* h_bc_in_sig_reg_2D_decays_[Counters::num_of_variable_pairs][Counters::num_of_processes][5][Counters::num_of_reco_channels];


};
#endif
