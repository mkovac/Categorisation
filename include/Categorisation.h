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
#include "TRandom3.h"

// Include classes
#include "Tree.h"
#include "Counters.h"
#include "Histograms.h"
#include "Utilities.h"
#include "bitops.h"
#include "cConstants.h"
#include "ROC.h"
#include "Category.h"
#include "Discriminants.h"


// Bools
#define APPLY_K_FACTORS 1
#define EXCLUDE_H2l2X 1
#define REQUIRE_EXACTLY_4_GOOD_LEPTONS      0
#define REQUIRE_AT_LEAST_5_GOOD_LEPTONS     0
#define REQUIRE_EXACTLY_5_GOOD_LEPTONS      0
#define REQUIRE_EXACTLY_6_GOOD_LEPTONS      0
#define REQUIRE_H_LEPTONS_ARE_IN_ETA_PT_ACC 0
#define REQUIRE_H_LEPTONS_ARE_GOOD          0

#define USEQGTAGGING 0

// Working points for MELA-only categorization
#define USEMASSDEPWPVBF2JMELA 0

// Working points for MELA+q/g categorization
#define USEMASSDEPWPVBF2JMELAQG 0


// MET
#define METCUT 100

// Baskets
#define BASKETLIST 12

// Jets
#define CSVv2M 0.8484
#define CSVv2L 0.5426
#define OFFICIALQGTAGGER 1
#define BTAGGINGSF 1 // 0 - none; 1 - central; 2 - up; 3 - down

// Checks
#define CHECKFROMCATEGORYCC 1


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
   void SumAssocAndHiggsLep();
   void FindGenChannel();
   void FindRecoChannel();
   void FillControlCounters();
   void ResetPerEventStuff();
   void SaveHistograms( TString );
   void DoLeptonMatching();
   void UseMatchingInfo();
   void FillHistograms();
   void MakeROCs();
   void Categorise();
   
   int FindCurrentProcess( TString );
   int FindCurrentAssocDecay();
   
   float CalculateFactorK( TString );
   
   bool IsSignal();


private:

   Histograms *histograms;
   ROC *roc;

   TFile *input_file_;
   TTree *input_tree_;

   TH1F* hist_counters_;
   
   int current_process_, current_final_state_, current_category_, current_assoc_dec_;
   float lumi_, k_factor_, partial_sample_weight_, m4l_min_, m4l_max_, delta_R_;
   double gen_sum_weights_, event_weight_, num_gen_events_;
   
   vector<TString> s_category_;
   
   int reco_ch_1, reco_ch_2, reco_ch_3;
   int gen_ch_1, gen_ch_2, gen_ch_3;
   
   // Working points
   float WP_VBF2j, WP_VBF1j, WP_WHh, WP_ZHh;
   
   // Matching statuses
   int H_lep_match_status   = -999;
   int all_lep_match_status = -999;
   int WH_lep_match_status  = -999;
   int ZH_lep_match_status  = -999;
   int ttH_lep_match_status = -999;
   int Z1_match_status      = -999;
   int Z2_match_status      = -999;
   
   // Current Decay
   int current_W_decay  = -999;
   int current_Z_decay  = -999;
   int current_tt_decay = -999;

   bool signal_region, pass_trigger, pass_trigger_no_1E;
   
   bool found_matching_ambiguity;
   
   int n_ones, n_ones_H_lep, n_ones_assoc_lep;
   
// Variables map
   map<Counters::variable, pair<float, bool>> variable_map;
   map<Counters::variable_pair, tuple<float, float, bool, bool> > variable_pair_map;
   
// Probabilities
   float Pvbf, Phjj;
 
// Discriminants
   float KD, D_2j_VBF_Hjj, D_1j_VBF_Hj, D_2j_WH_hadr_Hjj, D_2j_ZH_hadr_Hjj;

// Per event lepton variables
   vector<int>   gen_H_lep_id_;
   vector<float> gen_H_lep_pt_;
   vector<float> gen_H_lep_eta_;
   vector<float> gen_H_lep_phi_;
   vector<int>   gen_assoc_lep_id_;
   vector<float> gen_assoc_lep_pt_;
   vector<float> gen_assoc_lep_eta_;
   vector<float> gen_assoc_lep_phi_;
   
   vector<float> sorted_gen_H_lep_pt_;
   vector<float> sorted_gen_H_lep_abs_eta_;
   vector<float> sorted_cand_lep_pt_;
   vector<float> sorted_cand_lep_abs_eta_;

   
// Per event counters
   int n_gen_H_lep               = 0; // only e and mu
   int n_gen_H_lep_in_eta_acc    = 0;
   int n_gen_H_lep_in_pt_acc     = 0;
   int n_gen_H_lep_in_eta_pt_acc = 0;
   int n_gen_H_ele               = 0;
   int n_gen_H_mu                = 0;
   int n_gen_H_tau               = 0;
   int n_gen_H_LEP               = 0; // e, mu, and tau

   int n_gen_assoc_lep               = 0; // only e and mu
   int n_gen_assoc_lep_in_eta_acc    = 0;
   int n_gen_assoc_lep_in_pt_acc     = 0;
   int n_gen_assoc_lep_in_eta_pt_acc = 0;
   int n_gen_assoc_ele               = 0;
   int n_gen_assoc_mu                = 0;
   int n_gen_assoc_tau               = 0;
   int n_gen_assoc_LEP               = 0; // e, mu, and tau
   
   int n_gen_lep               = 0; // only e and mu
   int n_gen_lep_in_eta_acc    = 0;
   int n_gen_lep_in_pt_acc     = 0;
   int n_gen_lep_in_eta_pt_acc = 0;
   int n_gen_ele               = 0;
   int n_gen_mu                = 0;
   int n_gen_tau               = 0;
   int n_gen_LEP               = 0; // e, mu, and tau
   int n_gen_LEP_plus          = 0; // e, mu, and tau
   int n_gen_LEP_minus         = 0; // e, mu, and tau
   
   // Lepton matching counter
   map<TString, map<int, int>> counter_map; // move to pair???
   
   int n_reco_lep_matched_to_gen_H_lep[4]     = {0, 0, 0, 0};
   int n_cand_lep_matched_to_gen_H_lep[4]     = {0, 0, 0, 0};
   int n_reco_lep_matched_to_gen_assoc_lep[2] = {0, 0};
   int n_cand_lep_matched_to_gen_assoc_lep[2] = {0, 0};
   int n_gen_lep_matched_to_cand_lep[4]       = {0, 0, 0, 0};
   int n_gen_lep_matched_to_reco_lep[7]       = {0, 0, 0, 0, 0, 0, 0};
   int n_gen_H_lep_matched_to_Z1_lep[4]       = {0,0};
   int n_gen_H_lep_matched_to_Z2_lep[4]       = {0,0};

   
   bool  gen_H_lep_is_in_eta_acc[4];
   bool  gen_H_lep_is_in_pt_acc[4];
   bool  gen_H_lep_is_in_eta_pt_acc[4];
   bool  gen_assoc_lep_is_in_eta_acc[4];
   bool  gen_assoc_lep_is_in_pt_acc[4];
   bool  gen_assoc_lep_is_in_eta_pt_acc[4];
   
   
// Cuts impacting counters/histograms
   bool exactly_4_good_leptons_;
   bool at_least_5_good_leptons_;
   bool exactly_5_good_leptons_;
   bool exactly_6_good_leptons_;
   bool H_leptons_are_in_eta_pt_acc_;
   bool H_leptons_are_good_;
   
// Jets
   float jet_no_b_tag_;
   float jet_b_tag_;
   float vbf_2_jets;
   float vbf_lost_jet;
   
   float qg_is_default_;
   float qg_is_normal_;
   float jet_p_quark_[99];
   float jet_p_gluon_[99];
   float jet_p_g_over_p_q_[99];
   
   float jetPhi_[99];
   float jet_QG_likelihood_[99];
   float jet_QG_likelihood_raw_[99];
   
// Probabilities and discriminants
   float p_q_j1_p_q_j2;
   float p_g_j1_p_g_j2;
   
   float D_2j_qg ;
   float D_qg_j1_D_qg_j2;
   
   float p_quark;
   float p_gluon;
   float D_1j_qg;
   
   float D_2j_Mela_QG_VBF_Hjj;
   float D_2j_Mela_D_2j_QG_VBF_Hjj;
   
   float D_1j_Mela_QG_VBF_Hj;
   float D_1j_Mela_D_1j_QG_VBF_Hj;
   
   float D_2j_Mela_QG_WH_hadr_Hjj;
   float D_2j_Mela_D_2j_QG_WH_hadr_Hjj;
   
   float D_2j_Mela_QG_ZH_hadr_Hjj;
   float D_2j_Mela_D_2j_QG_ZH_hadr_Hjj;
   
   float D_2j_Mela_exp_QG_VBF_Hjj;
   
   float D_2j_Mela_sq_QG_VBF_Hjj;
   float D_2j_Mela_sqrt_QG_VBF_Hjj;
   float D_2j_Mela_cbrt_QG_VBF_Hjj;
   float D_2j_Mela_qrrt_QG_VBF_Hjj;
   float D_2j_Mela_qnrt_QG_VBF_Hjj;
   
   float D_1j_Mela_sqrt_QG_VBF_Hj;
   float D_1j_Mela_cbrt_QG_VBF_Hj;
   float D_1j_Mela_qrrt_QG_VBF_Hj;
   float D_1j_Mela_qnrt_QG_VBF_Hj;
   
   float D_2j_Mela_sqrt_QG_WH_hadr_Hjj;
   float D_2j_Mela_cbrt_QG_WH_hadr_Hjj ;
   float D_2j_Mela_qrrt_QG_WH_hadr_Hjj;
   float D_2j_Mela_qnrt_QG_WH_hadr_Hjj;
   
   float D_2j_Mela_sqrt_QG_ZH_hadr_Hjj;
   float D_2j_Mela_cbrt_QG_ZH_hadr_Hjj;
   float D_2j_Mela_qrrt_QG_ZH_hadr_Hjj;
   float D_2j_Mela_qnrt_QG_ZH_hadr_Hjj;
   
   
   
// Control counters
   int num_stored[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_gen_H_lep_in_eta_acc[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_gen_H_lep_in_eta_pt_acc[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_or_more_gen_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_or_more_reco_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_pass_triger[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_bc_in_sig_reg_and_pass_triger[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_bc_in_sig_reg_and_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_H_lep_in_eta_pt_acc_and_4_or_more_reco_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_H_lep_in_eta_pt_acc_and_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_H_lep_in_eta_pt_acc_and_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   int n_ev_4_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];

   float yield_stored[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_acc[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_4_or_more_gen_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_4_or_more_reco_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_bc_in_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_bc_in_sig_reg_and_pass_triger[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_bc_in_sig_reg_and_pass_triger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc_and_4_or_more_reco_lep[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc_and_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_H_lep_in_eta_pt_acc_pass_trig_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   float yield_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[Counters::num_of_processes][Counters::num_of_gen_ch];
   
// Counters
   int tot_n_bc_in_sig_reg_exactly_4_good_leps[Counters::num_of_processes][Counters::num_of_reco_ch];
   int tot_n_bc_in_sig_reg_at_least_5_good_leps[Counters::num_of_processes][Counters::num_of_reco_ch];
   int tot_n_bc_in_sig_reg_exactly_5_good_leps[Counters::num_of_processes][Counters::num_of_reco_ch];
   int tot_n_bc_in_sig_reg_exactly_6_good_leps[Counters::num_of_processes][Counters::num_of_reco_ch];
   int tot_n_bc_in_sig_reg_in_eta_pt_acc[Counters::num_of_processes][Counters::num_of_reco_ch];
   int tot_n_bc_in_sig_reg_H_leps_are_good[Counters::num_of_processes][Counters::num_of_reco_ch];

   int n_ev_bc_in_sr_all_4_leps_right_[Counters::num_of_processes][Counters::num_of_reco_ch];

   // WH
   int n_WH_events_[Counters::num_of_W_decays];
   int n_ev_H_leps_in_eta_pt_acc_WH_[Counters::num_of_W_decays];
   int n_ev_H_leps_are_good_WH_[Counters::num_of_W_decays];
   int n_ev_all_4_leps_right_WH_[Counters::num_of_W_decays];
   int n_Z1_daughters_from_HWH[Counters::num_of_Z1_match_statuses][Counters::num_of_W_decays];
   int n_Z2_daughters_from_HWH[Counters::num_of_Z2_match_statuses][Counters::num_of_W_decays];

   //ZH
   int n_ZH_events_[Counters::num_of_Z_decays];
   int n_ev_H_leps_in_eta_pt_acc_ZH_[Counters::num_of_Z_decays];
   int n_ev_H_leps_are_good_ZH_[Counters::num_of_Z_decays];
   int n_ev_all_4_leps_right_ZH_[Counters::num_of_Z_decays];
   int n_Z1_daughters_from_HZH[Counters::num_of_Z1_match_statuses][Counters::num_of_Z_decays];
   int n_Z2_daughters_from_HZH[Counters::num_of_Z2_match_statuses][Counters::num_of_Z_decays];
   
   //ttH
   int n_ttH_events_[Counters::num_of_Z_decays];
   int n_ev_H_leps_in_eta_pt_acc_ttH_[Counters::num_of_Z_decays];
   int n_ev_H_leps_are_good_ttH_[Counters::num_of_Z_decays];
   int n_ev_all_4_leps_right_ttH_[Counters::num_of_Z_decays];
   int n_Z1_daughters_from_HttH[Counters::num_of_Z1_match_statuses][Counters::num_of_Z_decays];
   int n_Z2_daughters_from_HttH[Counters::num_of_Z2_match_statuses][Counters::num_of_Z_decays];



   int n_ev_bc_in_sr_match_H_leps_[Counters::num_of_H_lep_match_statuses][Counters::num_of_processes][Counters::num_of_reco_ch];
   int n_ev_bc_in_sr_match_all_leps_[Counters::num_of_all_lep_match_statuses][Counters::num_of_processes][Counters::num_of_reco_ch];
   int n_ev_bc_in_sr_match_WH_leps_[Counters::num_of_WH_lep_match_statuses][Counters::num_of_processes][Counters::num_of_reco_ch];
   int n_ev_bc_in_sr_match_ZH_leps_[Counters::num_of_ZH_lep_match_statuses][Counters::num_of_processes][Counters::num_of_reco_ch];
   int n_ev_bc_in_sr_match_ttH_leps_[Counters::num_of_ttH_lep_match_statuses][Counters::num_of_processes][Counters::num_of_reco_ch];

   
// Reconstructed counters
   int num_of_events_with_bc[Counters::num_of_processes][Counters::num_of_reco_ch];
   int num_of_events_with_bc_in_sr[Counters::num_of_processes][Counters::num_of_reco_ch];
   float yield_of_events_with_bc[Counters::num_of_processes][Counters::num_of_reco_ch];
   float yield_of_events_with_bc_in_sr[Counters::num_of_processes][Counters::num_of_reco_ch];

   

// Higgs decay
 int gen_H_decay = -999;
   
};
#endif
