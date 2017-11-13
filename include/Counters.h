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


   static const int num_of_gen_channels        = MAX_NUM_OF_GEN_CHANNELS;
   static const int num_of_reco_channels       = MAX_NUM_OF_RECO_CHANNELS;
   static const int num_of_production_modes    = MAX_NUM_OF_PRODUCTION_MODES;
   static const int num_of_processes           = MAX_NUM_OF_PROCESSES;
   static const int num_of_final_states        = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_categories          = MAX_NUM_OF_CATEGORIES;
   static const int num_of_sorted_objects      = MAX_NUM_OF_SORTED_OBJECTS;
   static const int num_of_associated_decays   = MAX_NUM_OF_ASSOCIATED_DECAYS;

};
#endif
