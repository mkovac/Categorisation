#ifndef Counters_h
#define Counters_h

using namespace std;

class Counters
{

public:
	
	Counters();
	~Counters();
   
   enum final_state
   {
      fs4e    = 0,
      fs4mu   = 1,
      fs2e2mu = 2,
      fs2mu2e = 3,
      fs4l    = 4,
      MAX_NUM_OF_FINAL_STATES
   };
   
   enum generated_channel
   {
      gen_ch_4mu      = 0,
      gen_ch_4e       = 1,
      gen_ch_2e2mu    = 2,
      gen_ch_4tau     = 3,
      gen_ch_2e2tau   = 4,
      gen_ch_2mu2tau  = 5,
      gen_ch_other    = 6,
      gen_ch_e_mu     = 7,
      gen_ch_tau      = 8,
      gen_ch_2        = 9,
      gen_ch_3        = 10,
      MAX_NUM_OF_GEN_CHANNELS
   };
   
   enum reconstructed_channel
   {
      reco_ch_4mu      = 0,
      reco_ch_4e       = 1,
      reco_ch_2e2mu    = 2,
      reco_ch_all      = 3,
      reco_ch_1_def    = 4,
      reco_ch_2_def    = 5,
      reco_ch_3_def    = 6,
      MAX_NUM_OF_RECO_CHANNELS
   };
   
   enum category
   {
      untagged         = 0,
      VBF_1j_tagged    = 1,
      VBF_2j_tagged    = 2,
      VH_lepton_tagged = 3,
      VH_hadron_tagged = 4,
      ttH_tagged       = 5,
      VH_MET_tagged    = 6,
      inclusive        = 7,
      MAX_NUM_OF_CATEGORIES
   };
   
   enum production_mode
   {
      ggH = 0,
      VBF = 1,
      WH  = 2,
      ZH  = 3,
      ttH = 4,
      MAX_NUM_OF_PRODUCTION_MODES
   };
   
   enum process
   { 
      AllData = 0,
      H125    = 1,
      H125ggH = 2,
      H125VBF = 3,
      H125WH  = 4,
      H125ZH  = 5,
      H125bbH = 6,
      H125ttH = 7,
      qqZZ    = 8,
      ggZZ    = 9,
      ttbar   = 10,
      MAX_NUM_OF_PROCESSES
   };
   
   enum sort
   {
      first = 0,
      second = 1,
      third = 2,
      fourth  = 3,
      MAX_NUM_OF_SORTED_OBJECTS
   };

   enum associated_decay
   {
      WH_4_0  = 0,
      WH_4_1  = 1,
      ZH_4_0  = 2,
      ZH_4_2  = 3,
      ZH_2_2  = 4,
      bbH_0   = 5,
      bbH_1   = 6,
      ttH_4_0 = 7,
      ttH_4_1 = 8,
      ttH_4_2 = 9,
      ttH_2_2 = 10,
      MAX_NUM_OF_ASSOCIATED_DECAYS
   };

   enum variable
   {
      M4l,
      M4l2,
      MZ1,
      MZ2,
      Dkinbkg,
      DiJetFisher,
      Pvbf,
      Phjj,
      Pvbf1j,
      Phj,
      Pwhhadr,
      Pzhhadr,
      Pwhlept,
      Pzhlept,
      D2jVbfHjj,
      D1jVbfHj,
      D2jWHHadrHjj,
      D2jZHHadrHjj,
      Pqj1,
      Pgj1,
      Pqj1Pqj2,
      Pgj1Pgj2,
      D2jqg,
      Dqgj1Dqgj2,
      Pqj1VbfTopo,
      Pgj1VbfTopo,
      Pqj1Pqj2VbfTopo,
      Pgj1Pgj2VbfTopo,
      D2jqgVbfTopo,
      Pq,
      Pg,
      D1jqg,
      D2jMelaQGVbfHjj,
      D2jMelaD2jQGVbfHjj,
      D1jMelaQGVbfHj,
      D1jMelaD1jQGVbfHj,
      D2jMelaQGWHHadrHjj,
      D2jMelaD2jQGWHHadrHjj,
      D2jMelaQGZHHadrHjj,
      D2jMelaD2jQGZHHadrHjj,
      RatioPvbfPhjj,
      RatioPqj1Pqj2Pgj1Pgj2,
      RatioPvbfPqj1Pqj2PhjjPgj1Pgj2,
      D2jMelaExpQGVbfHjj,
      D2jMelaSqQGVbfHjj,
      D2jMelaSqrtQGVbfHjj,
      D2jMelaCbrtQGVbfHjj,
      D2jMelaQrrtQGVbfHjj,
      D2jMelaQnrtQGVbfHjj,
      D1jMelaSqrtQGVbfHj,
      D1jMelaCbrtQGVbfHj,
      D1jMelaQrrtQGVbfHj,
      D1jMelaQnrtQGVbfHj,
      D2jMelaSqrtQGWHHadrHjj,
      D2jMelaCbrtQGWHHadrHjj,
      D2jMelaQrrtQGWHHadrHjj,
      D2jMelaQnrtQGWHHadrHjj,
      D2jMelaSqrtQGZHHadrHjj,
      D2jMelaCbrtQGZHHadrHjj,
      D2jMelaQrrtQGZHHadrHjj,
      D2jMelaQnrtQGZHHadrHjj,
      Pt4l,
      NGenLep,
      NGenLepInEtaPtAcc,
      NGenLepNotInEtaPtAcc,
      NGenHLepNotInEtaPtAcc,
      NGenAssocLepNotInEtaPtAcc,
      NGenLepMinusNGoodLep,
      NGenLepInEtaPtAccMinusNGoodLep,
      NExtraLep,
      NExtraZ,
      NJets,
      NBtaggedJets,
      MET,
      MAX_NUMBER_OF_VARIABLES
   };
   
   enum variable_pair
   {
      M4l_vs_Dkinbkg,
      MZ2_vs_Dkinbkg,
      D2jVbfHjj_vs_D2jqg,
      MAX_NUMBER_OF_VARIABLE_PAIRS
   };
   
   enum H_lep_match_status
   {
      H_4         = 0,
      H_3         = 1,
      H_2         = 2,
      H_1         = 3,
      H_0         = 4,
      H_ambiguity = 5,
      MAX_NUMBER_OF_H_LEP_MATCH_STATUSES
   };
   
   // 40 = n_ones_H = 4, n_ones_assoc = 0
   enum all_lep_match_status
   {
      all_40        = 0,
      all_31        = 1,
      all_22        = 2,
      all_le_4      = 3,
      all_ambiguity = 4,
      MAX_NUMBER_OF_ALL_LEP_MATCH_STATUSES
   };
   
   // 40_40 = n_ones_H = 4, n_ones_assoc = 0, n_gen_H_lep = 4, n_gen_assoc_lep = 0
   enum WH_lep_match_status
   {
      WH_40_40     = 0,
      WH_40_41     = 1,
      WH_31_41     = 2,
      WH_le_4      = 3,
      WH_ambiguity = 4,
      MAX_NUMBER_OF_WH_LEP_MATCH_STATUSES
   };

   // 40_40 = n_ones_H = 4, n_ones_assoc = 0, n_gen_H_lep = 4, n_gen_assoc_lep = 0
   enum ZH_lep_match_status
   {
      ZH_40_40     = 0,
      ZH_40_42     = 1,
      ZH_31_42     = 2,
      ZH_22_42     = 3,
      ZH_22_22     = 4,
      ZH_le_4      = 5,
      ZH_ambiguity = 6,
      MAX_NUMBER_OF_ZH_LEP_MATCH_STATUSES
   };
   
   // 40_40 = n_ones_H = 4, n_ones_assoc = 0, n_gen_H_lep = 4, n_gen_assoc_lep = 0
   enum ttH_lep_match_status
   {
      ttH_40_40     = 0,
      ttH_40_41     = 1,
      ttH_31_41     = 2,
      ttH_40_42     = 3,
      ttH_31_42     = 3,
      ttH_22_42     = 5,
      ttH_22_22     = 6,
      ttH_le_4      = 7,
      ttH_ambiguity = 8,
      MAX_NUMBER_OF_TTH_LEP_MATCH_STATUSES
   };
   
   enum W_decays
   {
      W_dec_40 = 0,
      W_dec_41 = 1,
      MAX_NUMBER_OF_W_DECAYS
   };
   
   enum Z_decays
   {
      Z_dec_40 = 0,
      Z_dec_42 = 1,
      Z_dec_22 = 2,
      MAX_NUMBER_OF_Z_DECAYS
   };
   
   enum tt_decays
   {
      tt_dec_40 = 0,
      tt_dec_41 = 1,
      tt_dec_42 = 2,
      tt_dec_22 = 3,
      MAX_NUMBER_OF_tt_DECAYS
   };
   
   enum Z1_match_status
   {
      Z1_0 = 0,
      Z1_1 = 1,
      Z1_2 = 2,
      Z1_3 = 3,
      MAX_NUMBER_OF_Z1_MATCH_STATUSES
   };

   enum Z2_match_status
   {
      Z2_0 = 0,
      Z2_1 = 1,
      Z2_2 = 2,
      Z2_3 = 3,
      MAX_NUMBER_OF_Z2_MATCH_STATUSES
   };
   
   enum ROC
   {
      DiJetFisher_qqH_ggH,
      DiJetFisher_qqH_qqZZ,
      Pvbf_qqH_ggH,
      Pvbf_qqH_qqZZ,
      Phjj_qqH_ggH,
      Phjj_qqH_qqZZ,
      Phjj_WH_ggH,
      Phjj_ZH_ggH,
      Pvbf1j_qqH_ggH,
      Phj_qqH_ggH,
      Pwhhadr_WH_ggH,
      Pzhhadr_ZH_ggH,
      D2jVbfHjj_qqH_ggH,
      D2jVbfHjj_qqH_qqZZ,
      D1jVbfHj_qqH_ggH,
      D2jWHHadrHjj_WH_ggH,
      D2jZHHadrHjj_ZH_ggH,
      Pqj1Pqj2_qqH_ggH,
      Pgj1Pgj2_qqH_ggH,
      D2jqg_qqH_ggH,
      Dqgj1Dqgj2_qqH_ggH,
      Pqj1Pqj2_qqH_qqZZ,
      Pgj1Pgj2_qqH_qqZZ,
      D2jqg_qqH_qqZZ,
      Dqgj1Dqgj2_qqH_qqZZ,
      Pq_qqH_ggH,
      Pg_qqH_ggH,
      D1jqg_qqH_ggH,
      Pqj1Pqj2_WH_ggH,
      Pgj1Pgj2_WH_ggH,
      D2jqg_WH_ggH,
      Dqgj1Dqgj2_WH_ggH,
      Pqj1Pqj2_ZH_ggH,
      Pgj1Pgj2_ZH_ggH,
      D2jqg_ZH_ggH,
      Dqgj1Dqgj2_ZH_ggH,
      D2jMelaQGVbfHjj_qqH_ggH,
      D2jMelaD2jQGVbfHjj_qqH_ggH,
      D2jMelaQGVbfHjj_qqH_qqZZ,
      D2jMelaD2jQGVbfHjj_qqH_qqZZ,
      D1jMelaQGVbfHj_qqH_ggH,
      D1jMelaD1jQGVbfHj_qqH_ggH,
      D2jMelaQGWHHadrHjj_WH_ggH,
      D2jMelaD2jQGWHHadrHjj_WH_ggH,
      D2jMelaQGZHHadrHjj_ZH_ggH,
      D2jMelaD2jQGZHHadrHjj_ZH_ggH,
      RatioPvbfPhjj_qqH_ggH,
      RatioPqj1Pqj2Pgj1Pgj2_qqH_ggH,
      RatioPvbfPqj1Pqj2PhjjPgj1Pgj2_qqH_ggH,
      D2jMelaExpQGVbfHjj_qqH_ggH,
      D2jMelaSqQGVbfHjj_qqH_ggH,
      D2jMelaSqrtQGVbfHjj_qqH_ggH,
      D2jMelaCbrtQGVbfHjj_qqH_ggH,
      D2jMelaQrrtQGVbfHjj_qqH_ggH,
      D2jMelaQnrtQGVbfHjj_qqH_ggH,
      D2jMelaExpQGVbfHjj_qqH_qqZZ,
      D2jMelaSqQGVbfHjj_qqH_qqZZ,
      D2jMelaSqrtQGVbfHjj_qqH_qqZZ,
      D2jMelaCbrtQGVbfHjj_qqH_qqZZ,
      D2jMelaQrrtQGVbfHjj_qqH_qqZZ,
      D2jMelaQnrtQGVbfHjj_qqH_qqZZ,
      D1jMelaSqrtQGVbfHj_qqH_ggH,
      D1jMelaCbrtQGVbfHj_qqH_ggH,
      D1jMelaQrrtQGVbfHj_qqH_ggH,
      D1jMelaQnrtQGVbfHj_qqH_ggH,
      D2jMelaSqrtQGWHHadrHjj_WH_ggH,
      D2jMelaCbrtQGWHHadrHjj_WH_ggH,
      D2jMelaQrrtQGWHHadrHjj_WH_ggH,
      D2jMelaQnrtQGWHHadrHjj_WH_ggH,
      D2jMelaSqrtQGZHHadrHjj_ZH_ggH,
      D2jMelaCbrtQGZHHadrHjj_ZH_ggH,
      D2jMelaQrrtQGZHHadrHjj_ZH_ggH,
      D2jMelaQnrtQGZHHadrHjj_ZH_ggH,
      MAX_NUMBER_OF_ROCS
   };

   
   static const int num_of_gen_ch                 = MAX_NUM_OF_GEN_CHANNELS;
   static const int num_of_reco_ch                = MAX_NUM_OF_RECO_CHANNELS;
   static const int num_of_production_modes       = MAX_NUM_OF_PRODUCTION_MODES;
   static const int num_of_processes              = MAX_NUM_OF_PROCESSES;
   static const int num_of_final_states           = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_categories             = MAX_NUM_OF_CATEGORIES;
   static const int num_of_sorted_objects         = MAX_NUM_OF_SORTED_OBJECTS;
   static const int num_of_associated_decays      = MAX_NUM_OF_ASSOCIATED_DECAYS;
   static const int num_of_vars                   = MAX_NUMBER_OF_VARIABLES;
   static const int num_of_variable_pairs         = MAX_NUMBER_OF_VARIABLE_PAIRS;
   static const int num_of_H_lep_match_statuses   = MAX_NUMBER_OF_H_LEP_MATCH_STATUSES;
   static const int num_of_all_lep_match_statuses = MAX_NUMBER_OF_ALL_LEP_MATCH_STATUSES;
   static const int num_of_WH_lep_match_statuses  = MAX_NUMBER_OF_WH_LEP_MATCH_STATUSES;
   static const int num_of_ZH_lep_match_statuses  = MAX_NUMBER_OF_ZH_LEP_MATCH_STATUSES;
   static const int num_of_ttH_lep_match_statuses = MAX_NUMBER_OF_TTH_LEP_MATCH_STATUSES;
   static const int num_of_W_decays               = MAX_NUMBER_OF_W_DECAYS;
   static const int num_of_Z_decays               = MAX_NUMBER_OF_Z_DECAYS;
   static const int num_of_tt_decays              = MAX_NUMBER_OF_tt_DECAYS;
   static const int num_of_ROCs                   = MAX_NUMBER_OF_ROCS;

   static const int num_of_Z1_match_statuses      = MAX_NUMBER_OF_Z1_MATCH_STATUSES;
   static const int num_of_Z2_match_statuses      = MAX_NUMBER_OF_Z2_MATCH_STATUSES;   
};
#endif
