#ifndef Categorisation_h
#define Categorisation_h

// C++
#include <iostream>
#include <fstream>
#include <iomanip> // For setprecision
#include <vector>
#include <map>
#include <algorithm>

// ROOT
#include "TApplication.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TColor.h"

// Include classes
#include "Tree.h"
#include "Counters.h"
#include "Histograms.h"
#include "Utilities.h"
#include "bitops.h"


// BOOLS
#define APPLY_K_FACTORS 1
#define EXCLUDE_H2l2X 1


using namespace std;

class Categorisation: public Tree
{

public:
	
	Categorisation( float );
	~Categorisation();
   
   void MakeHistograms( TString );
   void CountAssocLep();
   void CountHiggsLep();
   void FillIdPtEtaPhi();
   void SumAssocAndH();
   void FindGenChannel();
   void FindRecoChannel();
   void FillControlCounters();
   void ResetPerEventStuff();
   void SaveHistograms( TString );
   void DoLeptonMatching();
   
   int FindCurrentProcess( TString );
   int FindCurrentAssocDecay();
   
   float CalculateFactorK( TString );
   
   bool IsSignal();


private:

   Histograms *histograms;

   TFile *p_input_file;
   TTree *p_input_tree;

   TH1F* p_hist_counters;
   
   int current_process_, current_final_state_, current_category_, current_assoc_dec_;
   float lumi_, k_factor_, partial_sample_weight_, m4l_min_, m4l_max_;
   double gen_sum_weights_, event_weight_, num_gen_events_;
   
   int reco_ch_1, reco_ch_2, reco_ch_3;
   int gen_ch_1, gen_ch_2, gen_ch_3;
   
   bool signal_region, pass_trigger, pass_trigger_no_1E;
   
// Per event lepton variables
   vector<int>   gen_H_lep_id_;
   vector<float> gen_H_lep_pt_;
   vector<float> gen_H_lep_eta_;
   vector<float> gen_H_lep_phi_;
   vector<int>   gen_assoc_lep_id_;
   vector<float> gen_assoc_lep_pt_;
   vector<float> gen_assoc_lep_eta_;
   vector<float> gen_assoc_lep_phi_;
   
// Lepton matching counters
   int n_reco_lep_matched_to_gen_H_lep[4]     = {0, 0, 0, 0};
   int n_cand_lep_matched_to_gen_H_lep[4]     = {0, 0, 0, 0};
   int n_reco_lep_matched_to_gen_assoc_lep[2] = {0, 0};
   int n_cand_lep_matched_to_gen_assoc_lep[2] = {0, 0};
   int n_gen_lep_matched_to_cand_lep[4]       = {0, 0, 0, 0};
   int n_gen_lep_matched_to_reco_lep[7]       = {0, 0, 0, 0, 0, 0, 0};
   int n_gen_H_lep_matched_to_Z1_lep[4]       = {0,0};
   int n_gen_H_lep_matched_to_Z2_lep[4]       = {0,0};

// Per event counters
   int n_gen_H_lep               = 0;
   int n_gen_H_lep_in_eta_acc    = 0;
   int n_gen_H_lep_in_pt_acc     = 0;
   int n_gen_H_lep_in_eta_pt_acc = 0;
   int n_gen_H_ele = 0;
   int n_gen_H_mu  = 0;
   int n_gen_H_tau = 0;
   int n_gen_H_LEP = 0; // including tau

   int n_gen_assoc_lep = 0;
   int n_gen_assoc_lep_in_eta_acc    = 0;
   int n_gen_assoc_lep_in_pt_acc     = 0;
   int n_gen_assoc_lep_in_eta_pt_acc = 0;
   int n_gen_assoc_ele = 0;
   int n_gen_assoc_mu  = 0;
   int n_gen_assoc_tau = 0;
   int n_gen_assoc_LEP = 0; // including tau
   
   int n_gen_LEP_plus  = 0; // including tau
   int n_gen_LEP_minus = 0; // including tau
   
   int n_gen_lep = 0;
   int n_gen_lep_in_eta_acc    = 0;
   int n_gen_lep_in_pt_acc     = 0;
   int n_gen_lep_in_eta_pt_acc = 0;
   int n_gen_ele = 0;
   int n_gen_mu  = 0;
   int n_gen_tau = 0;
   int n_gen_LEP = 0; // including tau
   
   bool  gen_H_lep_is_in_eta_acc[4];
   bool  gen_H_lep_is_in_pt_acc[4];
   bool  gen_H_lep_is_in_eta_pt_acc[4];
   bool  gen_assoc_lep_is_in_eta_acc[4];
   bool  gen_assoc_lep_is_in_pt_acc[4];
   bool  gen_assoc_lep_is_in_eta_pt_acc[4];
   
// Control counters
   int num_stored[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_acc[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_pt_acc[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_reco_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_pass_triger[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_bc_in_sig_reg_and_pass_triger[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_bc_in_sig_reg_and_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   int tot_n_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];

   float yield_stored[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_acc[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_pt_acc[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_reco_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_bc_in_sig_reg_and_pass_triger[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_bc_in_sig_reg_and_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   float yield_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[Counters::num_of_processes][Counters::num_of_gen_channels];
   
// Reconstructed counters
   int num_of_events_with_bc[Counters::num_of_processes][Counters::num_of_reco_channels];
   int num_of_events_with_bc_sr[Counters::num_of_processes][Counters::num_of_reco_channels];
   float yield_of_events_with_bc[Counters::num_of_processes][Counters::num_of_reco_channels];
   float yield_of_events_with_bc_sr[Counters::num_of_processes][Counters::num_of_reco_channels];

   
   vector<float> sorted_gen_H_lep_pt_, sorted_gen_H_lep_abs_eta_, sorted_cand_lep_pt_, sorted_cand_lep_abs_eta_;

// Higgs decay
 int gen_H_decay = -999;
   
};
#endif
