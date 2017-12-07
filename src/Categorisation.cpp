// Include classes
#include "Categorisation.h"

// Constructor
//============================================================
Categorisation::Categorisation( float lumi ):Tree()
{
   lumi_ = lumi;
   current_process_     = -999;
   current_final_state_ = -999;
   current_category_    = -999;
   k_factor_ = 1;
   
   // Jets
   jet_b_tag_    = 0;
   jet_no_b_tag_ = 0;
   vbf_2_jets    = 0;
   vbf_lost_jet  = 0;

   
   m4l_min_ = 118.;
   m4l_max_ = 130.;
   
   gen_ch_1 = -999;
   
   reco_ch_1 = Counters::reco_ch_1_def;
   reco_ch_2 = Counters::reco_ch_2_def;
   reco_ch_3 = Counters::reco_ch_3_def;
   
   histograms = new Histograms(lumi_);
   
   signal_region = (ZZsel >= 90);
   pass_trigger = test_bit_16(trigWord, 0);
   pass_trigger_no_1E = test_bit_16(trigWord, 8);
   

}
//============================================================



// Destructor
//====================
Categorisation::~Categorisation()
{
}
//====================



//=============================================================
void Categorisation::MakeHistograms( TString input_file_name )
{

   p_input_file = TFile::Open(input_file_name);
   
   p_hist_counters = (TH1F*)p_input_file->Get("ZZTree/Counters");
   num_gen_events_ = (double)p_hist_counters->GetBinContent(1);
   gen_sum_weights_ = (double)p_hist_counters->GetBinContent(40);
   
   p_input_tree = (TTree*)p_input_file->Get("ZZTree/candTree");
   Init( p_input_tree, input_file_name );
   
   // Find current process
   current_process_ = FindCurrentProcess(input_file_name);
   
   if ( fChain == 0 ) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      // Reset per-event counters
      ResetPerEventStuff();
   
      // K factors
      if ( APPLY_K_FACTORS ) k_factor_ = CalculateFactorK(input_file_name);
   
      // Final event weight
      event_weight_ = (lumi_ * 1000 * xsec * k_factor_ * overallEventWeight)/gen_sum_weights_;
   
//   cout << event_weight_ << endl;
   


      FillIdPtEtaPhi();
      CountHiggsLep();
      CountAssocLep();
      SumAssocAndH();
      FindGenChannel();
      FindRecoChannel();
      FillControlCounters();
      
      if ( EXCLUDE_H2l2X )
      {
         if ( n_gen_H_lep != 4 && IsSignal() ) continue;
      }
      
      // Find associated decay
      current_assoc_dec_ = FindCurrentAssocDecay();
      
      
      // Successive selection steps
      num_of_events_with_bc[current_process_][reco_ch_1]++;
      num_of_events_with_bc[current_process_][reco_ch_2]++;
      yield_of_events_with_bc[current_process_][reco_ch_1] += event_weight_;
      yield_of_events_with_bc[current_process_][reco_ch_2] += event_weight_;

      if ( !signal_region ) continue;

      if ( ZZMass < m4l_min_ || ZZMass > m4l_max_ ) continue;

      num_of_events_with_bc_sr[current_process_][reco_ch_1]++;
      num_of_events_with_bc_sr[current_process_][reco_ch_2]++;
      yield_of_events_with_bc_sr[current_process_][reco_ch_1] += event_weight_;
      yield_of_events_with_bc_sr[current_process_][reco_ch_2] += event_weight_;
      
      // Match
      DoLeptonMatching();
      
      // Cuts impacting counters and histograms
      exactly_4_good_leptons_      = nExtraLep == 0;
      at_least_5_good_leptons_     = nExtraLep >= 1;
      exactly_5_good_leptons_      = nExtraLep == 1;
      exactly_6_good_leptons_      = nExtraLep == 2;
      H_leptons_are_in_eta_pt_acc_ = n_gen_H_lep_in_eta_pt_acc == 4;
      H_leptons_are_good_          = !found_matching_ambiguity && n_ones == 4;

      if ( REQUIRE_EXACTLY_4_GOOD_LEPTONS && !exactly_4_good_leptons_ )           continue;
      if ( REQUIRE_AT_LEAST_5_GOOD_LEPTONS && !at_least_5_good_leptons_ )         continue;
      if ( REQUIRE_EXACTLY_5_GOOD_LEPTONS && !exactly_5_good_leptons_ )           continue;
      if ( REQUIRE_EXACTLY_6_GOOD_LEPTONS && !exactly_6_good_leptons_ )           continue;
      if ( REQUIRE_H_LEPTONS_ARE_IN_ETA_PT_ACC && !H_leptons_are_in_eta_pt_acc_ ) continue;
      if ( REQUIRE_H_LEPTONS_ARE_GOOD && !H_leptons_are_good_ )                   continue;

      FillHistograms();

   } // end events loop
}
//=====================================================



//=================================================================
int Categorisation::FindCurrentProcess( TString input_file_name )
{
   
   int current_process = -999;
   
   // Assign dataset to correct process
   if ( input_file_name.Contains("AllData") )        current_process = Counters::AllData;
   if ( input_file_name.Contains("ggH125") )         current_process = Counters::H125ggH;
   if ( input_file_name.Contains("VBFH125") )        current_process = Counters::H125VBF;
   if ( input_file_name.Contains("WplusH125") )      current_process = Counters::H125WH;
   if ( input_file_name.Contains("WminusH125") )     current_process = Counters::H125WH;
   if ( input_file_name.Contains("ZH125") )          current_process = Counters::H125ZH;
   if ( input_file_name.Contains("bbH125") )         current_process = Counters::H125bbH;
   if ( input_file_name.Contains("ttH125") )         current_process = Counters::H125ttH;
   if ( input_file_name.Contains("ZZTo4l") )         current_process = Counters::qqZZ;
   if ( input_file_name.Contains("ggTo4e") )         current_process = Counters::ggZZ;
   if ( input_file_name.Contains("ggTo4mu") )        current_process = Counters::ggZZ;
   if ( input_file_name.Contains("ggTo4tau") )       current_process = Counters::ggZZ;
   if ( input_file_name.Contains("ggTo2e2mu") )      current_process = Counters::ggZZ;
   if ( input_file_name.Contains("ggTo2e2tau") )     current_process = Counters::ggZZ;
   if ( input_file_name.Contains("ggTo2mu2tau") )    current_process = Counters::ggZZ;
   // End assign dataset to correct process
   
   return current_process;
}
//=================================================================



//=================================================================
bool Categorisation::IsSignal()
{
   
   bool is_signal = 0;
   
   // Assign dataset to correct process
   if ( current_process_ == Counters::AllData )   is_signal = 0;
   if ( current_process_ == Counters::H125ggH )   is_signal = 1;
   if ( current_process_ == Counters::H125VBF )   is_signal = 1;
   if ( current_process_ == Counters::H125WH )    is_signal = 1;
   if ( current_process_ == Counters::H125WH)     is_signal = 1;
   if ( current_process_ == Counters::H125ZH )    is_signal = 1;
   if ( current_process_ == Counters::H125bbH )   is_signal = 1;
   if ( current_process_ == Counters::H125ttH )   is_signal = 1;
   if ( current_process_ == Counters::qqZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   if ( current_process_ == Counters::ggZZ )      is_signal = 0;
   // End assign dataset to correct process
   
   return is_signal;
}
//=================================================================



//=========================================
int Categorisation::FindCurrentAssocDecay()
{
   
   int current_assoc_dec = -999;
   
   if ( current_process_ == Counters::H125WH )
   {
      if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )     current_assoc_dec = Counters::WH_4_0;
      else if (n_gen_H_lep == 4 && n_gen_assoc_lep == 1 ) current_assoc_dec = Counters::WH_4_1;
   }
   
   if ( current_process_ == Counters::H125ZH )
   {
      if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )      current_assoc_dec = Counters::ZH_4_0;
      else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 2 ) current_assoc_dec = Counters::ZH_4_2;
      else if ( n_gen_H_lep == 2 && n_gen_assoc_lep == 2 ) current_assoc_dec = Counters::ZH_2_2;
   }
   
   if ( current_process_ == Counters::H125ttH )
   {
      if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )      current_assoc_dec = Counters::ttH_4_0;
      else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 1 ) current_assoc_dec = Counters::ttH_4_1;
      else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 2 ) current_assoc_dec = Counters::ttH_4_2;
      else if ( n_gen_H_lep == 2 && n_gen_assoc_lep == 2 ) current_assoc_dec = Counters::ttH_2_2;
   }
   
   bool trigger_error = (current_process_ == Counters::H125WH || current_process_ == Counters::H125ZH || current_process_ == Counters::H125ttH)
                        && current_assoc_dec == -999;
   if ( trigger_error )
   {
      cout << "[ERROR] Cannot find current associated decay!" << endl;
   }
   
   return current_assoc_dec;
}
//=========================================



//=================================================================
float Categorisation::CalculateFactorK( TString input_file_name )
{

   float k_factor = 1;
   
   if ( input_file_name.Contains("ZZTo4l"))
   {
      k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
   }
   else if ( input_file_name.Contains("ggTo"))
   {
      k_factor = KFactor_QCD_ggZZ_Nominal;
   }
   return k_factor;
}
//=================================================================


//====================================
void Categorisation::FillIdPtEtaPhi()
{
   gen_H_lep_id_.push_back(GenLep1Id);
   gen_H_lep_id_.push_back(GenLep2Id);
   gen_H_lep_id_.push_back(GenLep3Id);
   gen_H_lep_id_.push_back(GenLep4Id);

   gen_H_lep_pt_.push_back(GenLep1Pt);
   gen_H_lep_pt_.push_back(GenLep2Pt);
   gen_H_lep_pt_.push_back(GenLep3Pt);
   gen_H_lep_pt_.push_back(GenLep4Pt);

   gen_H_lep_eta_.push_back(GenLep1Eta);
   gen_H_lep_eta_.push_back(GenLep2Eta);
   gen_H_lep_eta_.push_back(GenLep3Eta);
   gen_H_lep_eta_.push_back(GenLep4Eta);

   gen_H_lep_phi_.push_back(GenLep1Phi);
   gen_H_lep_phi_.push_back(GenLep2Phi);
   gen_H_lep_phi_.push_back(GenLep3Phi);
   gen_H_lep_phi_.push_back(GenLep4Phi);

   gen_assoc_lep_id_.push_back(GenAssocLep1Id);
   gen_assoc_lep_id_.push_back(GenAssocLep2Id);

   gen_assoc_lep_pt_.push_back(GenAssocLep1Pt);
   gen_assoc_lep_pt_.push_back(GenAssocLep2Pt);
   
   gen_assoc_lep_eta_.push_back(GenAssocLep1Eta);
   gen_assoc_lep_eta_.push_back(GenAssocLep2Eta);

   gen_assoc_lep_phi_.push_back(GenAssocLep1Phi);
   gen_assoc_lep_phi_.push_back(GenAssocLep2Phi);
}
//====================================


//===================================
void Categorisation::CountHiggsLep()
{
   for ( int i_gen_H_lep = 0; i_gen_H_lep < 4; i_gen_H_lep++ )
   {
      if ( abs(gen_H_lep_id_.at(i_gen_H_lep)) == 11 )
      {
         n_gen_H_ele++;
         n_gen_H_lep++;
         n_gen_H_LEP++;
        
        if ( gen_H_lep_id_.at(i_gen_H_lep) > 0 )
        {
           n_gen_LEP_plus++;
        }
        else
        {
           n_gen_LEP_minus++;
        }
        
        gen_H_lep_is_in_eta_acc[i_gen_H_lep]    = fabs(gen_H_lep_eta_.at(i_gen_H_lep)) < 2.5;
        gen_H_lep_is_in_pt_acc[i_gen_H_lep]     = gen_H_lep_pt_.at(i_gen_H_lep) > 7;
        gen_H_lep_is_in_eta_pt_acc[i_gen_H_lep] = gen_H_lep_is_in_eta_acc[i_gen_H_lep] && gen_H_lep_is_in_pt_acc[i_gen_H_lep];

        if ( gen_H_lep_is_in_eta_acc[i_gen_H_lep] )    n_gen_H_lep_in_eta_acc++;
        if ( gen_H_lep_is_in_pt_acc[i_gen_H_lep] )     n_gen_H_lep_in_pt_acc++;
        if ( gen_H_lep_is_in_eta_pt_acc[i_gen_H_lep] ) n_gen_H_lep_in_eta_pt_acc++;
      }
      else if ( abs(gen_H_lep_id_.at(i_gen_H_lep)) == 13 )
      {
         n_gen_H_mu++;
         n_gen_H_lep++;
         n_gen_H_LEP++;
         
         if ( gen_H_lep_id_.at(i_gen_H_lep) > 0 )
         {
            n_gen_LEP_plus++;
         }
         else
         {
            n_gen_LEP_minus++;
         }
         
         gen_H_lep_is_in_eta_acc[i_gen_H_lep]    = fabs(gen_H_lep_eta_.at(i_gen_H_lep)) < 2.2;
         gen_H_lep_is_in_pt_acc[i_gen_H_lep]     = gen_H_lep_pt_.at(i_gen_H_lep) > 5;
         gen_H_lep_is_in_eta_pt_acc[i_gen_H_lep] = gen_H_lep_is_in_eta_acc[i_gen_H_lep] && gen_H_lep_is_in_pt_acc[i_gen_H_lep];
         
         if ( gen_H_lep_is_in_eta_acc[i_gen_H_lep] )    n_gen_H_lep_in_eta_acc++;
         if ( gen_H_lep_is_in_pt_acc[i_gen_H_lep] )     n_gen_H_lep_in_pt_acc++;
         if ( gen_H_lep_is_in_eta_pt_acc[i_gen_H_lep] ) n_gen_H_lep_in_eta_pt_acc++;
      }
      else if ( abs(gen_H_lep_id_.at(i_gen_H_lep)) == 15 )
      {
         n_gen_H_tau++;
         n_gen_H_LEP++;
         
         if ( gen_H_lep_id_.at(i_gen_H_lep) > 0 )
         {
            n_gen_LEP_plus++;
         }
         else
         {
            n_gen_LEP_minus++;
         }
      }
   } // end for
}
//===================================


//===================================
void Categorisation::CountAssocLep()
{
   for ( int i_gen_assoc_lep = 0; i_gen_assoc_lep < 2; i_gen_assoc_lep++ )
   {
      if ( abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 11 )
      {
         n_gen_assoc_ele++;
         n_gen_assoc_lep++;
         n_gen_assoc_LEP++;
        
        if ( gen_assoc_lep_id_.at(i_gen_assoc_lep) > 0 )
        {
           n_gen_LEP_plus++;
        }
        else
        {
           n_gen_LEP_minus++;
        }
        
        gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep]    = fabs(gen_assoc_lep_eta_.at(i_gen_assoc_lep)) < 2.5;
        gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep]     = gen_assoc_lep_pt_.at(i_gen_assoc_lep) > 7;
        gen_assoc_lep_is_in_eta_pt_acc[i_gen_assoc_lep] = gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep] && gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep];

        if ( gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep] )    n_gen_assoc_lep_in_eta_acc++;
        if ( gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep] )     n_gen_assoc_lep_in_pt_acc++;
        if ( gen_assoc_lep_is_in_eta_pt_acc[i_gen_assoc_lep] ) n_gen_assoc_lep_in_eta_pt_acc++;
      }
      else if ( abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 13 )
      {
         n_gen_assoc_mu++;
         n_gen_assoc_lep++;
         n_gen_assoc_LEP++;
         
         if ( gen_assoc_lep_id_.at(i_gen_assoc_lep) > 0 )
         {
            n_gen_LEP_plus++;
         }
         else
         {
            n_gen_LEP_minus++;
         }
         
         gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep]    = fabs(gen_assoc_lep_eta_.at(i_gen_assoc_lep)) < 2.2;
         gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep]     = gen_assoc_lep_pt_.at(i_gen_assoc_lep) > 5;
         gen_assoc_lep_is_in_eta_pt_acc[i_gen_assoc_lep] = gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep] && gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep];
         
         if ( gen_assoc_lep_is_in_eta_acc[i_gen_assoc_lep] )    n_gen_assoc_lep_in_eta_acc++;
         if ( gen_assoc_lep_is_in_pt_acc[i_gen_assoc_lep] )     n_gen_assoc_lep_in_pt_acc++;
         if ( gen_assoc_lep_is_in_eta_pt_acc[i_gen_assoc_lep] ) n_gen_assoc_lep_in_eta_pt_acc++;
      }
      else if ( abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 15 )
      {
         n_gen_assoc_tau++;
         n_gen_assoc_LEP++;
         
         if ( gen_assoc_lep_id_.at(i_gen_assoc_lep) > 0 )
         {
            n_gen_LEP_plus++;
         }
         else
         {
            n_gen_LEP_minus++;
         }
      }
   } // end for
}
//===================================


//=================================
void Categorisation::SumAssocAndH()
{
   n_gen_lep               = n_gen_H_lep + n_gen_assoc_lep;
   n_gen_lep_in_eta_acc    = n_gen_H_lep_in_eta_acc + n_gen_assoc_lep_in_eta_acc;
   n_gen_lep_in_pt_acc     = n_gen_H_lep_in_pt_acc + n_gen_assoc_lep_in_pt_acc;
   n_gen_lep_in_eta_pt_acc = n_gen_H_lep_in_eta_pt_acc + n_gen_assoc_lep_in_eta_pt_acc;
   n_gen_ele               = n_gen_H_ele + n_gen_assoc_ele;
   n_gen_mu                = n_gen_H_mu + n_gen_assoc_mu;
   n_gen_tau               = n_gen_H_tau + n_gen_assoc_tau;
   n_gen_LEP               = n_gen_H_LEP + n_gen_assoc_LEP;
}
//=================================


//===================================
void Categorisation::FindGenChannel()
{
   int temp_gen_ch = -999;
   
   gen_ch_2 = 9; // channels 0, 1, or 2
   gen_ch_3 = 10; // channels 3, 4, or 5

   if      ( n_gen_H_mu  == 4 )                     temp_gen_ch = Counters::gen_ch_4mu;
   else if ( n_gen_H_ele == 4 )                     temp_gen_ch = Counters::gen_ch_4e;
   else if ( n_gen_H_ele == 2 && n_gen_H_mu == 2 )  temp_gen_ch = Counters::gen_ch_2e2mu;
   else if ( n_gen_H_tau == 4 )                     temp_gen_ch = Counters::gen_ch_4tau;
   else if ( n_gen_H_ele == 2 && n_gen_H_tau == 2 ) temp_gen_ch = Counters::gen_ch_2e2tau;
   else if ( n_gen_H_mu  == 2 && n_gen_H_tau == 2 ) temp_gen_ch = Counters::gen_ch_2mu2tau;
   else                                             temp_gen_ch = Counters::gen_ch_other;
   
   gen_ch_1 = temp_gen_ch;
   
   if ( 0 <= temp_gen_ch && temp_gen_ch <= 2 ) gen_ch_2 = 7;
   if ( 3 <= temp_gen_ch && temp_gen_ch <= 5 ) gen_ch_2 = 8;
}
//===================================


//====================================
void Categorisation::FindRecoChannel()
{

   int n_cand_ele = -999;
   int n_cand_mu  = -999;

   for ( int i_cand_lep = 0; i_cand_lep < 4; i_cand_lep++ )
   {
      if ( abs(LepLepId->at(i_cand_lep)) == 11 ) n_cand_ele++;
      if ( abs(LepLepId->at(i_cand_lep)) == 13 ) n_cand_mu++;
   }
      if      ( n_cand_ele == 0 && n_cand_mu == 4 ) reco_ch_1 = Counters::reco_ch_4mu;
      else if ( n_cand_ele == 4 && n_cand_mu == 0 ) reco_ch_1 = Counters::reco_ch_4e;
      else if ( n_cand_ele == 2 && n_cand_mu == 2 ) reco_ch_1 = Counters::reco_ch_2e2mu;
   
      bool all_reco_ch = reco_ch_1 == Counters::reco_ch_4mu || reco_ch_1 == Counters::reco_ch_4e || reco_ch_1 == Counters::reco_ch_2e2mu;
      if ( all_reco_ch ) reco_ch_2 = Counters::reco_ch_all;
}
//====================================


//========================================
void Categorisation::FillControlCounters()
{
   num_stored[current_process_][gen_ch_1]++;
   num_stored[current_process_][gen_ch_2]++;
   num_stored[current_process_][gen_ch_3]++;
   
   yield_stored[current_process_][gen_ch_1] += event_weight_;
   yield_stored[current_process_][gen_ch_2] += event_weight_;
   yield_stored[current_process_][gen_ch_3] += event_weight_;

//   cout << n_gen_H_lep_in_eta_acc << endl;


   if ( n_gen_H_lep_in_eta_acc == 4 )
   {
      tot_n_gen_H_lep_in_eta_acc[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_acc[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_acc[current_process_][gen_ch_3]++;
   
      yield_gen_H_lep_in_eta_acc[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_acc[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_acc[current_process_][gen_ch_3] += event_weight_;
      
      sorted_gen_H_lep_pt_ = gen_H_lep_pt_;
      sort(sorted_gen_H_lep_pt_.begin(), sorted_gen_H_lep_pt_.end(), std::greater<>());
      
      for ( int i_sort = 0; i_sort < 4; i_sort++ )
      {
//         cout << sorted_gen_H_lep_pt_.at(i_sort) << endl;
         histograms->FillPt(sorted_gen_H_lep_pt_.at(i_sort), event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3, i_sort);
      }
   } // end if
   
   
   if ( n_gen_H_lep_in_pt_acc == 4 )
   {
      for ( vector<float>::iterator it = gen_H_lep_eta_.begin(); it != gen_H_lep_eta_.end(); it++ )
      {
         sorted_gen_H_lep_abs_eta_.push_back(abs(*it));
      } // end for
      
      sort(sorted_gen_H_lep_abs_eta_.begin(), sorted_gen_H_lep_abs_eta_.end(), std::greater<>());

      for ( int i_sort = 0; i_sort < 4; i_sort++ )
      {
//         cout << sorted_gen_H_lep_abs_eta_.at(i_sort) << endl;
         histograms->FillAbsEta(sorted_gen_H_lep_abs_eta_.at(i_sort), event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3, i_sort);
      } // end for

   } // end if



   if ( n_gen_lep_in_eta_pt_acc == 4 )
   {
      tot_n_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_3]++;
      
      yield_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc[current_process_][gen_ch_3] += event_weight_;
      
//      cout << NRecoMu << " " << gen_ch_1 << endl;
      
      histograms->FillRecoEleN(NRecoEle, event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3);
      histograms->FillRecoMuN(NRecoMu, event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3);
      histograms->FillRecoLepN(NRecoMu + NRecoEle, event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3);
   }
   
   if ( n_gen_lep >= 4 )
   {
      tot_n_gen_lep[current_process_][gen_ch_1]++;
      tot_n_gen_lep[current_process_][gen_ch_2]++;
      tot_n_gen_lep[current_process_][gen_ch_3]++;

      yield_gen_lep[current_process_][gen_ch_1] += event_weight_;
      yield_gen_lep[current_process_][gen_ch_2] += event_weight_;
      yield_gen_lep[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( NRecoMu + NRecoEle >= 4 )
   {
      tot_n_reco_lep[current_process_][gen_ch_1]++;
      tot_n_reco_lep[current_process_][gen_ch_2]++;
      tot_n_reco_lep[current_process_][gen_ch_3]++;

      yield_reco_lep[current_process_][gen_ch_1] += event_weight_;
      yield_reco_lep[current_process_][gen_ch_2] += event_weight_;
      yield_reco_lep[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( signal_region )
   {
      tot_n_bc_in_sig_reg[current_process_][gen_ch_1]++;
      tot_n_bc_in_sig_reg[current_process_][gen_ch_2]++;
      tot_n_bc_in_sig_reg[current_process_][gen_ch_3]++;

      yield_bc_in_sig_reg[current_process_][gen_ch_1] += event_weight_;
      yield_bc_in_sig_reg[current_process_][gen_ch_2] += event_weight_;
      yield_bc_in_sig_reg[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( pass_trigger )
   {
      tot_n_pass_triger[current_process_][gen_ch_1]++;
      tot_n_pass_triger[current_process_][gen_ch_2]++;
      tot_n_pass_triger[current_process_][gen_ch_3]++;

      yield_pass_trigger[current_process_][gen_ch_1] += event_weight_;
      yield_pass_trigger[current_process_][gen_ch_2] += event_weight_;
      yield_pass_trigger[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( pass_trigger_no_1E )
   {
      tot_n_pass_triger_no_1E[current_process_][gen_ch_1]++;
      tot_n_pass_triger_no_1E[current_process_][gen_ch_2]++;
      tot_n_pass_triger_no_1E[current_process_][gen_ch_3]++;

      yield_pass_triger_no_1E[current_process_][gen_ch_1] += event_weight_;
      yield_pass_triger_no_1E[current_process_][gen_ch_2] += event_weight_;
      yield_pass_triger_no_1E[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( signal_region && pass_trigger )
   {
      tot_n_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_1]++;
      tot_n_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_2]++;
      tot_n_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_3]++;

      yield_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_1] += event_weight_;
      yield_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_2] += event_weight_;
      yield_bc_in_sig_reg_and_pass_triger[current_process_][gen_ch_3] += event_weight_;
   
      sorted_cand_lep_pt_ = *LepPt;
      sort(sorted_cand_lep_pt_.begin(), sorted_cand_lep_pt_.end(), std::greater<>());
      
      for ( int i_sort = 0; i_sort < 4; i_sort++ )
      {
//         cout << sorted_cand_lep_pt_.at(i_sort) << endl;
         histograms->FillPtReco(sorted_cand_lep_pt_.at(i_sort), event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3, i_sort);
      }

      for ( vector<float>::iterator it = LepEta->begin(); it != LepEta->end(); it++ )
      {
         sorted_cand_lep_abs_eta_.push_back(abs(*it));
      } // end for
      
      sort(sorted_cand_lep_abs_eta_.begin(), sorted_cand_lep_abs_eta_.end(), std::greater<>());

      for ( int i_sort = 0; i_sort < 4; i_sort++ )
      {
//         cout << sorted_cand_lep_abs_eta_.at(i_sort) << endl;
         histograms->FillAbsEtaReco(sorted_cand_lep_abs_eta_.at(i_sort), event_weight_, current_process_, gen_ch_1, gen_ch_2, gen_ch_3, i_sort);
      } // end for
   
   } // end if


   if ( signal_region && pass_trigger_no_1E )
   {
      tot_n_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_1]++;
      tot_n_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_2]++;
      tot_n_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_3]++;

      yield_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_1] += event_weight_;
      yield_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_2] += event_weight_;
      yield_bc_in_sig_reg_and_pass_triger_no_1E[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( n_gen_H_lep_in_eta_pt_acc == 4 && (NRecoMu + NRecoEle) >= 4 )
   {
      tot_n_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_3]++;

      yield_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_4_reco_lep[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( n_gen_H_lep_in_eta_pt_acc == 4 && signal_region )
   {
      tot_n_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_3]++;

      yield_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_bc_in_sig_reg[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( n_gen_H_lep_in_eta_pt_acc == 4 && pass_trigger )
   {
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_3]++;

      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger[current_process_][gen_ch_3] += event_weight_;
   }
   
   
   if ( n_gen_H_lep_in_eta_pt_acc == 4 && pass_trigger_no_1E )
   {
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_1]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_2]++;
      tot_n_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_3]++;

      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_1] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_2] += event_weight_;
      yield_gen_H_lep_in_eta_pt_acc_and_pass_trigger_no_1E[current_process_][gen_ch_3] += event_weight_;
   }

   if ( n_gen_H_lep_in_eta_pt_acc == 4 && signal_region && pass_trigger )
   {
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_1]++;
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_2]++;
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_3]++;

      yield_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_1] += event_weight_;
      yield_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_2] += event_weight_;
      yield_gen_lep_in_eta_pt_acc_pass_trig_sig_reg[current_process_][gen_ch_3] += event_weight_;
   }
   
   if ( n_gen_H_lep_in_eta_pt_acc == 4 && signal_region && pass_trigger_no_1E )
   {
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_1]++;
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_2]++;
      tot_n_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_3]++;

      yield_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_1] += event_weight_;
      yield_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_2] += event_weight_;
      yield_gen_lep_in_eta_pt_acc_pass_trig_no_1E_sig_reg[current_process_][gen_ch_3] += event_weight_;
   }
   
   histograms->FillPtEtaH(GenHPt, GenHRapidity, event_weight_, current_process_);
}
//========================================



//=====================================
void Categorisation::DoLeptonMatching()
{
   for ( int i_gen_H_lep = 0; i_gen_H_lep < 4; i_gen_H_lep++ )
   {
      if ( abs(gen_H_lep_id_.at(i_gen_H_lep)) == 11 || abs(gen_H_lep_id_.at(i_gen_H_lep)) == 13 )
      {
         for ( int i_cand_lep = 0; i_cand_lep < 4; i_cand_lep++ )
         {
            delta_R_ = Utilities::DeltaR(gen_H_lep_eta_.at(i_gen_H_lep), gen_H_lep_phi_.at(i_gen_H_lep),
                                         LepEta->at(i_cand_lep), LepPhi->at(i_cand_lep));
            if ( delta_R_ < 0.1 )
            {
               counter_map["n_reco_lep_matched_to_gen_H_lep"][i_gen_H_lep]++;
               counter_map["n_cand_lep_matched_to_gen_H_lep"][i_gen_H_lep]++;
               counter_map["n_gen_lep_matched_to_reco_lep"][i_cand_lep]++;
               counter_map["n_gen_lep_matched_to_cand_lep"][i_cand_lep]++;

//               n_reco_lep_matched_to_gen_H_lep[i_gen_H_lep]++;
//               n_cand_lep_matched_to_gen_H_lep[i_gen_H_lep]++;
//               n_gen_lep_matched_to_reco_lep[i_cand_lep]++;
//               n_gen_lep_matched_to_cand_lep[i_cand_lep]++;
         
               if ( i_cand_lep < 2 )
               {
                  counter_map["n_gen_H_lep_matched_to_Z1_lep"][i_cand_lep]++;
//                  n_gen_H_lep_matched_to_Z1_lep[i_cand_lep]++;
               }
               else
               {
                  counter_map["n_gen_H_lep_matched_to_Z2_lep"][i_cand_lep - 2]++;
//                  n_gen_H_lep_matched_to_Z2_lep[i_cand_lep - 2]++;
               } // end else
            } // end if
         } // end for
         
         for ( int i_extra_lep = 0; i_extra_lep < nExtraLep; i_extra_lep++ )
         {
            delta_R_ = Utilities::DeltaR(gen_H_lep_eta_.at(i_gen_H_lep), gen_H_lep_phi_.at(i_gen_H_lep),
                                         ExtraLepEta->at(i_extra_lep), ExtraLepPhi->at(i_extra_lep));
            if ( delta_R_ < 0.1 )
            {
               counter_map["n_reco_lep_matched_to_gen_H_lep"][i_gen_H_lep]++;
               counter_map["n_gen_lep_matched_to_reco_lep"][4 + i_extra_lep]++;

//               n_reco_lep_matched_to_gen_H_lep[i_gen_H_lep]++;
//               n_gen_lep_matched_to_reco_lep[4 + i_extra_lep]++;

            } // end if
         } // end for
      } // end if
   } // end for
   
   for ( int i_gen_assoc_lep = 0; i_gen_assoc_lep < 2; i_gen_assoc_lep++ )
   {
      if ( abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 11 || abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 13 )
      {
         for ( int i_cand_lep = 0; i_cand_lep < 4; i_cand_lep++ )
         {
            delta_R_ = Utilities::DeltaR(gen_H_lep_eta_.at(i_gen_assoc_lep), gen_H_lep_phi_.at(i_gen_assoc_lep),
                                         LepEta->at(i_cand_lep), LepPhi->at(i_cand_lep));
            if ( delta_R_ < 0.1 )
            {
               counter_map["n_reco_lep_matched_to_gen_assoc_lep"][i_gen_assoc_lep]++;
               counter_map["n_cand_lep_matched_to_gen_assoc_lep"][i_gen_assoc_lep]++;
               counter_map["n_gen_lep_matched_to_reco_lep"][i_cand_lep]++;
               counter_map["n_gen_lep_matched_to_cand_lep"][i_cand_lep]++;
               
//               n_reco_lep_matched_to_gen_assoc_lep[i_gen_assoc_lep]++;
//               n_cand_lep_matched_to_gen_assoc_lep[i_gen_assoc_lep]++;
//               n_gen_lep_matched_to_reco_lep[i_cand_lep]++;
//               n_gen_lep_matched_to_cand_lep[i_cand_lep]++;

            } // end if
         } // end for
         
         for ( int i_extra_lep = 0; i_extra_lep < nExtraLep; i_extra_lep++ )
         {
            delta_R_ = Utilities::DeltaR(gen_assoc_lep_eta_.at(i_gen_assoc_lep), gen_assoc_lep_phi_.at(i_gen_assoc_lep),
                                         ExtraLepEta->at(i_extra_lep), ExtraLepPhi->at(i_extra_lep));
            if ( delta_R_ < 0.1 )
            {
               counter_map["n_reco_lep_matched_to_gen_assoc_lep"][i_gen_assoc_lep]++;
               counter_map["n_gen_lep_matched_to_reco_lep"][4 + i_extra_lep]++;
               
//               n_reco_lep_matched_to_gen_assoc_lep[i_gen_assoc_lep]++;
//               n_gen_lep_matched_to_reco_lep[4 + i_extra_lep]++;
               
            } // end if
         } // end for
      } // end if
   } // end for
}
//=====================================


//====================================
void Categorisation::UseMatchingInfo()
{
   found_matching_ambiguity = false;
   
   n_ones = 0;
   n_ones_H_lep = 0;
   n_ones_assoc_lep = 0;
   
   
   for ( int i_gen_H_lep = 0; i_gen_H_lep < 4; i_gen_H_lep++ )
   {
      if ( counter_map["n_reco_lep_matched_to_gen_H_lep"][i_gen_H_lep] > 1 )
      {
         found_matching_ambiguity = true;
         break;
      }
   }
   
   for ( int i_gen_assoc_lep = 0; i_gen_assoc_lep < 2; i_gen_assoc_lep++ )
   {
      if ( counter_map["n_reco_lep_matched_to_gen_assoc_lep"][i_gen_assoc_lep] > 1 )
      {
         found_matching_ambiguity = true;
         break;
      }
   }
   
   for ( int i_reco_lep = 0; i_reco_lep < 4 + nExtraLep; i_reco_lep++ )
   {
      if ( counter_map["n_gen_lep_matched_to_reco_lep"][i_reco_lep] > 1 )
      {
         found_matching_ambiguity = true;
         break;
      }
   }




   for ( int i_gen_H_lep = 0; i_gen_H_lep < 4; i_gen_H_lep++ )
   {
      if ( counter_map["n_reco_lep_matched_to_gen_H_lep"][i_gen_H_lep] == 1 ) n_ones++;
   }
   
   for ( int i_gen_H_lep = 0; i_gen_H_lep < 4; i_gen_H_lep++ )
   {
      if ( counter_map["n_cand_lep_matched_to_gen_H_lep"][i_gen_H_lep] == 1 ) n_ones_H_lep++;
   }
    
   for ( int i_gen_assoc_lep = 0; i_gen_assoc_lep < 2; i_gen_assoc_lep++ )
   {
      if ( counter_map["n_cand_lep_matched_to_gen_assoc_lep"][i_gen_assoc_lep] == 1 ) n_ones_assoc_lep++;
   }
   
   
   int current_match_H_lep_status = -1;
   int current_match_all_lep_status = -1;
   
   if ( found_matching_ambiguity)
   {
      current_match_H_lep_status = 5;
   }
   else
   {
      if ( n_ones_H_lep == 4 ) current_match_H_lep_status = 0;
      if ( n_ones_H_lep == 3 ) current_match_H_lep_status = 1;
      if ( n_ones_H_lep == 2 ) current_match_H_lep_status = 2;
      if ( n_ones_H_lep == 1 ) current_match_H_lep_status = 3;
      if ( n_ones_H_lep == 0 ) current_match_H_lep_status = 4;
   }
   
   if ( found_matching_ambiguity )
   {
      current_match_all_lep_status = 4;
   }
   else
   {
      if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 ) current_match_all_lep_status = 0;
      if ( n_ones_H_lep == 3 && n_ones_assoc_lep == 1 ) current_match_all_lep_status = 1;
      if ( n_ones_H_lep == 2 && n_ones_assoc_lep == 2 ) current_match_all_lep_status = 2;
      if ( n_ones_H_lep + n_ones_assoc_lep < 4 )        current_match_all_lep_status = 3;
   }
   
   // WH, ZH and ttH match status
   int current_match_WH_status = -1;
   int current_match_ZH_status = -1;
   int current_match_ttH_status = -1;
   
   // WH
   if ( current_process_ == Counters::H125WH )
   {
      if ( found_matching_ambiguity )
      {
         current_match_WH_status = 4;
      }
      else
      {
         if ( n_ones_H_lep + n_ones_assoc_lep < 4 )
         {
            current_match_WH_status = 3;
         }
         else
         {
            if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_WH_status = 0;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 1 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_WH_status = 1;
               }
               else if ( n_ones_H_lep == 3 && n_ones_assoc_lep == 1 )
               {
                  current_match_WH_status = 2;
               }
               else cout << "[ERROR} n_ones" << endl;
            }
            else cout << "[ERROR] n_gen" << endl;
         }
      }
   }
   
   // ZH
   if ( current_process_ == Counters::H125ZH )
   {
      if ( found_matching_ambiguity )
      {
         current_match_ZH_status = 6;
      }
      else
      {
         if ( n_ones_H_lep + n_ones_assoc_lep < 4 )
         {
            current_match_ZH_status = 5;
         }
         else
         {
            if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_ZH_status = 0;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 2 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_ZH_status = 1;
               }
               else if ( n_ones_H_lep == 3 && n_ones_assoc_lep == 1 )
               {
                  current_match_ZH_status = 2;
               }
               else if ( n_ones_H_lep == 2 && n_ones_assoc_lep == 2 )
               {
                  current_match_ZH_status = 3;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 2 && n_gen_assoc_lep == 2 )
            {
               if ( n_ones_H_lep == 2 && n_ones_assoc_lep == 2 )
               {
                  current_match_ZH_status = 4;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else cout << "[ERROR] n_gen" << endl;
         }
      }
   }
   
   // ttH
   if ( current_process_ == Counters::H125ttH )
   {
      if ( found_matching_ambiguity )
      {
         current_match_ttH_status = 8;
      }
      else
      {
         if ( n_ones_H_lep + n_ones_assoc_lep < 4 )
         {
            current_match_ttH_status = 7;
         }
         else
         {
            if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 0 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_ttH_status = 0;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 1 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_ttH_status = 1;
               }
               else if ( n_ones_H_lep == 3 && n_ones_assoc_lep == 1 )
               {
                  current_match_ttH_status = 2;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 4 && n_gen_assoc_lep == 2 )
            {
               if ( n_ones_H_lep == 4 && n_ones_assoc_lep == 0 )
               {
                  current_match_ttH_status = 3;
               }
               else if ( n_ones_H_lep == 3 && n_ones_assoc_lep == 1 )
               {
                  current_match_ttH_status = 4;
               }
               else if ( n_ones_H_lep == 2 && n_ones_assoc_lep == 2 )
               {
                  current_match_ttH_status = 5;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else if ( n_gen_H_lep == 2 && n_gen_assoc_lep == 2 )
            {
               if ( n_ones_H_lep == 2 && n_ones_assoc_lep == 2 )
               {
                  current_match_ttH_status = 6;
               }
               else cout << "[ERROR] n_ones" << endl;
            }
            else cout << "[ERROR] n_gen" << endl;
         }
      }
   }

   int current_Z1_match_status = -1;
   int current_Z2_match_status = -1;
   
   if ( n_gen_H_lep_matched_to_Z1_lep[0] > 1 || n_gen_H_lep_matched_to_Z1_lep[1] > 1 )
   {
      current_Z1_match_status = 3;
   }
   else
   {
      current_Z1_match_status = 2 - (n_gen_H_lep_matched_to_Z1_lep[0] + n_gen_H_lep_matched_to_Z1_lep[1]);
   }
   
   if ( n_gen_H_lep_matched_to_Z2_lep[0] > 1 || n_gen_H_lep_matched_to_Z2_lep[1] > 1 )
   {
      current_Z2_match_status = 3;
   }
   else
   {
      current_Z2_match_status = 2 - (n_gen_H_lep_matched_to_Z2_lep[0] + n_gen_H_lep_matched_to_Z2_lep[1]);
   }
}
//====================================



//====================================
void Categorisation::FillHistograms()
{

   if ( exactly_4_good_leptons_ )
   {
      tot_n_bc_in_sig_reg_exactly_4_good_leps[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_exactly_4_good_leps[current_process_][reco_ch_2]++;
   }
   
   if ( at_least_5_good_leptons_ )
   {
      tot_n_bc_in_sig_reg_at_least_5_good_leps[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_at_least_5_good_leps[current_process_][reco_ch_2]++;
   }
   
   if ( exactly_5_good_leptons_ )
   {
      tot_n_bc_in_sig_reg_exactly_5_good_leps[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_exactly_5_good_leps[current_process_][reco_ch_2]++;
   }
   
      if ( exactly_6_good_leptons_ )
   {
      tot_n_bc_in_sig_reg_exactly_6_good_leps[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_exactly_6_good_leps[current_process_][reco_ch_2]++;
   }
   
   if ( H_leptons_are_in_eta_pt_acc_ )
   {
      tot_n_bc_in_sig_reg_in_eta_pt_acc[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_in_eta_pt_acc[current_process_][reco_ch_2]++;
   }
   
   if ( H_leptons_are_good_ )
   {
      tot_n_bc_in_sig_reg_H_leps_are_good[current_process_][reco_ch_1]++;
      tot_n_bc_in_sig_reg_H_leps_are_good[current_process_][reco_ch_2]++;
   }
   
   
   // Jets
   float c_QG_unoff = 1./1000.;
   int n_jets_B_tagged_redone = 0;
   int n_jets_B_tagged_loose_redone = 0;
   
   for ( int i_jet = 0; i_jet < nCleanedJetsPt30; i_jet++ )
   {
      if ( JetBTagger->at(i_jet) > CSVv2M ) n_jets_B_tagged_redone++;
      if ( JetBTagger->at(i_jet) > CSVv2L ) n_jets_B_tagged_loose_redone++;
   
      if ( i_jet < 2 )
      {
         if ( JetIsBtagged->at(i_jet) )
         {
            jet_b_tag_ += event_weight_;
         }
         else jet_no_b_tag_ += event_weight_;
      }
      
      
      jet_QG_likelihood_[i_jet] = JetQGLikelihood->at(i_jet);
   
      if ( JetQGLikelihood->at(i_jet) < 0. && i_jet < 2 )
      {
         TRandom3 rand;
         rand.SetSeed(abs(static_cast<int>(sin(JetPhi->at(i_jet))*100000)));
         jet_QG_likelihood_[i_jet] = rand.Uniform();
         
         qg_is_default_ += event_weight_;
      }
      else
      {
         qg_is_normal_ += event_weight_;
      }
   
      jet_p_quark_[i_jet] = 0.; //JetPQuark->at(j);//1000000. * JetPQuark->at(j);
      jet_p_gluon_[i_jet] = 0.; //JetPGluon->at(j);//1000. * JetPGluon->at(j);
      jet_p_g_over_p_q_[i_jet] = OFFICIALQGTAGGER ? (1./jet_QG_likelihood_[i_jet] - 1.) : (c_QG_unoff*jet_p_gluon_[i_jet]/jet_p_quark_[i_jet]);
   }
   
   
//   Fix this line
//   if ( !BTAGGINGSF && nCleanedJetsPt30BTagged!= n_jets_B_tagged_redone ) cout<<"ERROR : inconsistency in number of b-tagged jets"<<endl;
   
   if ( nCleanedJets >= 2 )
   {
      vbf_2_jets += event_weight_;
   }
   else vbf_lost_jet += event_weight_;
   
   float vbf_lost_jet = 0.; // Why?
   float c_VBF_2j = getDVBF2jetsConstant(ZZMass);
   float c_VBF_1j = getDVBF1jetConstant(ZZMass);
   float c_WH = getDWHhConstant(ZZMass);
   float c_ZH = getDZHhConstant(ZZMass);
   
   float c_WH_lept = 100000000000.;
   float c_ZH_lept = 100000000.;
   float p_WH_lept = p_LepWH_SIG_ghw1_1_JHUGen/c_WH_lept;
   float p_ZH_lept = p_LepZH_SIG_ghz1_1_JHUGen/c_ZH_lept;
   
   float KD = 1/(1 + getDbkgkinConstant(Z1Flav*Z2Flav, ZZMass)*p_QQB_BKG_MCFM/p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
   
   float D_2j_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
   float D_1j_VBF_Hj  = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal));
   
   float D_2j_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
   float D_2j_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
   
   
   if ( nCleanedJets >= 2 )
   {
      p_q_j1_p_q_j2 = 1000*1000*jet_p_quark_[0]*jet_p_quark_[1]/(c_QG_unoff*c_QG_unoff);
      p_g_j1_p_g_j2 = 1000*1000*jet_p_gluon_[0]*jet_p_gluon_[1];
   
      D_2j_qg = 1/(1 + jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1]);
      D_qg_j1_D_qg_j2 = 1/(1 + jet_p_g_over_p_q_[0]) * 1/(1 + jet_p_g_over_p_q_[1]);
   
      D_2j_Mela_QG_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1]);
   
      D_2j_Mela_QG_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal*jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1]);
      D_2j_Mela_QG_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal*jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1]);

      D_2j_Mela_exp_QG_VBF_Hjj = 1./(1.+ c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Exp(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1]));
   
      D_2j_Mela_sq_QG_VBF_Hjj   = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Power(jet_p_g_over_p_q_[0] * jet_p_g_over_p_q_[1],2));
      D_2j_Mela_sqrt_QG_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Power(jet_p_g_over_p_q_[0] * jet_p_g_over_p_q_[1],1/2));
      D_2j_Mela_cbrt_QG_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Power(jet_p_g_over_p_q_[0] * jet_p_g_over_p_q_[1],1/3));
      D_2j_Mela_qrrt_QG_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Power(jet_p_g_over_p_q_[0] * jet_p_g_over_p_q_[1],1/4));
      D_2j_Mela_qnrt_QG_VBF_Hjj = 1/(1 + c_VBF_2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal*TMath::Power(jet_p_g_over_p_q_[0] * jet_p_g_over_p_q_[1],1/5));
   
      D_2j_Mela_sqrt_QG_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/2));
      D_2j_Mela_cbrt_QG_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/3));
      D_2j_Mela_qrrt_QG_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/4));
      D_2j_Mela_qnrt_QG_WH_hadr_Hjj = 1/(1 + c_WH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/5));
   
      D_2j_Mela_sqrt_QG_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/2));
      D_2j_Mela_cbrt_QG_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/3));
      D_2j_Mela_qrrt_QG_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/4));
      D_2j_Mela_qnrt_QG_ZH_hadr_Hjj = 1/(1 + c_ZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1],1/5));
   }
   else
   {
      p_q_j1_p_q_j2 = -2;
      p_g_j1_p_g_j2 = -2;
   
      D_2j_qg = -2;
      D_qg_j1_D_qg_j2 = -2;
   
      D_2j_Mela_QG_VBF_Hjj = -2;
   
      D_2j_Mela_QG_WH_hadr_Hjj = -2;
      D_2j_Mela_QG_ZH_hadr_Hjj = -2;

      D_2j_Mela_exp_QG_VBF_Hjj = -2;
   
      D_2j_Mela_sq_QG_VBF_Hjj   = -2;
      D_2j_Mela_sqrt_QG_VBF_Hjj = -2;
      D_2j_Mela_cbrt_QG_VBF_Hjj = -2;
      D_2j_Mela_qrrt_QG_VBF_Hjj = -2;
      D_2j_Mela_qnrt_QG_VBF_Hjj = -2;
   
      D_2j_Mela_sqrt_QG_WH_hadr_Hjj = -2;
      D_2j_Mela_cbrt_QG_WH_hadr_Hjj = -2;
      D_2j_Mela_qrrt_QG_WH_hadr_Hjj = -2;
      D_2j_Mela_qnrt_QG_WH_hadr_Hjj = -2;
   
      D_2j_Mela_sqrt_QG_ZH_hadr_Hjj = -2;
      D_2j_Mela_cbrt_QG_ZH_hadr_Hjj = -2;
      D_2j_Mela_qrrt_QG_ZH_hadr_Hjj = -2;
      D_2j_Mela_qnrt_QG_ZH_hadr_Hjj = -2;
   }
   
   if ( nCleanedJets >= 1 )
   {
      p_quark = 1000 * jet_p_quark_[0]/c_QG_unoff;
      p_gluon = 1000 * jet_p_gluon_[0];
      D_1j_qg = 1/(1 + jet_p_g_over_p_q_[0]);
      
      D_1j_Mela_QG_VBF_Hj = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal)*jet_p_g_over_p_q_[0]);
   
      D_1j_Mela_sqrt_QG_VBF_Hj = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jet_p_g_over_p_q_[0],1/2));
      D_1j_Mela_cbrt_QG_VBF_Hj = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jet_p_g_over_p_q_[0],1/3));
      D_1j_Mela_qrrt_QG_VBF_Hj = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jet_p_g_over_p_q_[0],1/4));
      D_1j_Mela_qnrt_QG_VBF_Hj = 1/(1 + (c_VBF_1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jet_p_g_over_p_q_[0],1/5));
   }
   else
   {
      p_quark = -2;
      p_gluon = -2;
      D_1j_qg = -2;
      
      D_1j_Mela_QG_VBF_Hj = -2;
   
      D_1j_Mela_sqrt_QG_VBF_Hj = -2;
      D_1j_Mela_cbrt_QG_VBF_Hj = -2;
      D_1j_Mela_qrrt_QG_VBF_Hj = -2;
      D_1j_Mela_qnrt_QG_VBF_Hj = -2;
   }
   
   D_2j_Mela_D_2j_QG_VBF_Hjj = D_2j_VBF_Hjj * D_2j_qg;
   D_1j_Mela_D_1j_QG_VBF_Hj  = D_1j_VBF_Hj * D_1j_qg;
   
   D_2j_Mela_D_2j_QG_WH_hadr_Hjj = D_2j_WH_hadr_Hjj * D_2j_qg;
   D_2j_Mela_D_2j_QG_ZH_hadr_Hjj = D_2j_ZH_hadr_Hjj * D_2j_qg;
   
   
   cout << "[INFO] Filling variable map..." << endl;
   
   variables_map[Counters::M4l].first  = ZZMass;
   variables_map[Counters::M4l].second = 1;

   variables_map[Counters::M4l2].first  = ZZMass;
   variables_map[Counters::M4l2].second = 1;
   
   variables_map[Counters::MZ1].first  = Z1Mass;
   variables_map[Counters::MZ1].second = 1;
   
   variables_map[Counters::MZ2].first  = Z2Mass;
   variables_map[Counters::MZ2].second = 1;
   
   variables_map[Counters::Dkinbkg].first  = KD;
   variables_map[Counters::Dkinbkg].second = 1;
   
   variables_map[Counters::DiJetFisher].first = DiJetFisher;
   variables_map[Counters::DiJetFisher].second = (nCleanedJets >= 2);

   variables_map[Counters::Pvbf].first  = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal/c_VBF_2j;
   variables_map[Counters::Pvbf].second = (nCleanedJets >= 2);
   
   variables_map[Counters::Phjj].first  = p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
   variables_map[Counters::Phjj].second = (nCleanedJets >= 2);
   
   variables_map[Counters::Pvbf1j].first  = p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal/c_VBF_1j;
   variables_map[Counters::Pvbf1j].second = nCleanedJets == 1;

   variables_map[Counters::Phj].first  = p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
   variables_map[Counters::Phj].second = nCleanedJets == 1;
   
   variables_map[Counters::Pwhhadr].first  = p_HadWH_SIG_ghw1_1_JHUGen_JECNominal/c_WH;
   variables_map[Counters::Pwhhadr].second = (nCleanedJets >= 2);
   
   variables_map[Counters::Pzhhadr].first  = p_HadZH_SIG_ghz1_1_JHUGen_JECNominal/c_ZH;
   variables_map[Counters::Pzhhadr].second = (nCleanedJets >= 2);
   
   variables_map[Counters::Pwhlept].first  = p_WH_lept;
   variables_map[Counters::Pwhlept].second = nExtraLep >= 1;

   variables_map[Counters::Pzhlept].first  = p_ZH_lept;
   variables_map[Counters::Pzhlept].second = nExtraZ >= 1;
   
   variables_map[Counters::D2jVbfHjj].first  = D_2j_VBF_Hjj;
   variables_map[Counters::D2jVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D1jVbfHj].first  = D_1j_VBF_Hj;
   variables_map[Counters::D1jVbfHj].second = nCleanedJets == 1;

   variables_map[Counters::D2jWHHadrHjj].first  = D_2j_WH_hadr_Hjj;
   variables_map[Counters::D2jWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jZHHadrHjj].first  = D_2j_ZH_hadr_Hjj;
   variables_map[Counters::D2jZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::Pqj1].first  = p_quark;
   variables_map[Counters::Pqj1].second = (nCleanedJets >= 2);

   variables_map[Counters::Pgj1].first  = p_gluon;
   variables_map[Counters::Pgj1].second = (nCleanedJets >= 2);

   variables_map[Counters::Pqj1Pqj2].first  = p_q_j1_p_q_j2;
   variables_map[Counters::Pqj1Pqj2].second = (nCleanedJets >= 2);

   variables_map[Counters::Pgj1Pgj2].first  = p_g_j1_p_g_j2;
   variables_map[Counters::Pgj1Pgj2].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jqg].first  = D_2j_qg;
   variables_map[Counters::D2jqg].second = (nCleanedJets >= 2);

   variables_map[Counters::Dqgj1Dqgj2].first  = D_qg_j1_D_qg_j2;
   variables_map[Counters::Dqgj1Dqgj2].second = (nCleanedJets >= 2);

   variables_map[Counters::Pqj1VbfTopo].first  = p_quark;
   variables_map[Counters::Pqj1VbfTopo].second = (nCleanedJets >= 2 && D_2j_VBF_Hjj > 0.5);

   variables_map[Counters::Pgj1VbfTopo].first  = p_gluon;
   variables_map[Counters::Pgj1VbfTopo].second = (nCleanedJets >= 2 && D_2j_VBF_Hjj > 0.5);

   variables_map[Counters::Pqj1Pqj2VbfTopo].first  = p_q_j1_p_q_j2;
   variables_map[Counters::Pqj1Pqj2VbfTopo].second = (nCleanedJets >= 2 && D_2j_VBF_Hjj > 0.5);

   variables_map[Counters::Pgj1Pgj2VbfTopo].first  = p_g_j1_p_g_j2;
   variables_map[Counters::Pgj1Pgj2VbfTopo].second = (nCleanedJets >= 2 && D_2j_VBF_Hjj > 0.5);

   variables_map[Counters::D2jqgVbfTopo].first  = D_2j_qg;
   variables_map[Counters::D2jqgVbfTopo].second = (nCleanedJets >= 2 && D_2j_VBF_Hjj > 0.5);

   variables_map[Counters::Pq].first  = p_quark;
   variables_map[Counters::Pq].second = (nCleanedJets == 1);

   variables_map[Counters::Pg].first  = p_gluon;
   variables_map[Counters::Pg].second = (nCleanedJets == 1);

   variables_map[Counters::D1jqg].first  = D_1j_qg;
   variables_map[Counters::D1jqg].second = (nCleanedJets == 1);

   variables_map[Counters::D2jMelaQGVbfHjj].first  = D_2j_Mela_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaQGVbfHjj].second = (nCleanedJets >= 2);
   
   variables_map[Counters::D2jMelaD2jQGVbfHjj].first  = D_2j_Mela_D_2j_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaD2jQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D1jMelaQGVbfHj].first  = D_1j_Mela_QG_VBF_Hj;
   variables_map[Counters::D1jMelaQGVbfHj].second = (nCleanedJets == 1);

   variables_map[Counters::D1jMelaD1jQGVbfHj].first  = D_1j_Mela_D_1j_QG_VBF_Hj;
   variables_map[Counters::D1jMelaD1jQGVbfHj].second = (nCleanedJets == 1);

   variables_map[Counters::D2jMelaQGWHHadrHjj].first  = D_2j_Mela_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaQGWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaD2jQGWHHadrHjj].first  = D_2j_Mela_D_2j_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaD2jQGWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQGZHHadrHjj].first  = D_2j_Mela_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaD2jQGZHHadrHjj].first  = D_2j_Mela_D_2j_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaD2jQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::RatioPvbfPhjj].first  = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal/p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
   variables_map[Counters::RatioPvbfPhjj].second = (nCleanedJets >= 2);

   variables_map[Counters::RatioPqj1Pqj2Pgj1Pgj2].first  = jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1];
   variables_map[Counters::RatioPqj1Pqj2Pgj1Pgj2].second = (nCleanedJets >= 2);

   variables_map[Counters::RatioPvbfPqj1Pqj2PhjjPgj1Pgj2].first  = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal/p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal*jet_p_g_over_p_q_[0]*jet_p_g_over_p_q_[1];
   variables_map[Counters::RatioPvbfPqj1Pqj2PhjjPgj1Pgj2].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaExpQGVbfHjj].first  = D_2j_Mela_exp_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaExpQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaSqQGVbfHjj].first  = D_2j_Mela_sq_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaSqQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaSqrtQGVbfHjj].first  = D_2j_Mela_sqrt_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaSqrtQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaCbrtQGVbfHjj].first  = D_2j_Mela_cbrt_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaCbrtQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQrrtQGVbfHjj].first  = D_2j_Mela_qrrt_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaQrrtQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQnrtQGVbfHjj].first  = D_2j_Mela_qnrt_QG_VBF_Hjj;
   variables_map[Counters::D2jMelaQnrtQGVbfHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D1jMelaSqrtQGVbfHj].first  = D_1j_Mela_sqrt_QG_VBF_Hj;
   variables_map[Counters::D1jMelaSqrtQGVbfHj].second = (nCleanedJets == 1);
   
   variables_map[Counters::D1jMelaCbrtQGVbfHj].first  = D_1j_Mela_cbrt_QG_VBF_Hj;
   variables_map[Counters::D1jMelaCbrtQGVbfHj].second = (nCleanedJets == 1);
   
   variables_map[Counters::D1jMelaQrrtQGVbfHj].first  = D_1j_Mela_qrrt_QG_VBF_Hj;
   variables_map[Counters::D1jMelaQrrtQGVbfHj].second = (nCleanedJets == 1);

   variables_map[Counters::D1jMelaQnrtQGVbfHj].first  = D_1j_Mela_qnrt_QG_VBF_Hj;
   variables_map[Counters::D1jMelaQnrtQGVbfHj].second = (nCleanedJets == 1);

   variables_map[Counters::D2jMelaSqrtQGWHHadrHjj].first  = D_2j_Mela_sqrt_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaSqrtQGWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaCbrtQGWHHadrHjj].first  = D_2j_Mela_cbrt_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaCbrtQGWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQrrtQGWHHadrHjj].first  = D_2j_Mela_qrrt_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaQrrtQGWHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQnrtQGWHHadrHjj].first  = D_2j_Mela_qnrt_QG_WH_hadr_Hjj;
   variables_map[Counters::D2jMelaQnrtQGWHHadrHjj].second = (nCleanedJets >= 2);
   
   variables_map[Counters::D2jMelaSqrtQGZHHadrHjj].first  = D_2j_Mela_sqrt_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaSqrtQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaCbrtQGZHHadrHjj].first  = D_2j_Mela_cbrt_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaCbrtQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQrrtQGZHHadrHjj].first  = D_2j_Mela_qrrt_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaQrrtQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::D2jMelaQnrtQGZHHadrHjj].first  = D_2j_Mela_qnrt_QG_ZH_hadr_Hjj;
   variables_map[Counters::D2jMelaQnrtQGZHHadrHjj].second = (nCleanedJets >= 2);

   variables_map[Counters::Pt4l].first  = ZZPt;
   variables_map[Counters::Pt4l].second = 1;

   variables_map[Counters::NGenLep].first  = (float)n_gen_lep;
   variables_map[Counters::NGenLep].second = 1;

   variables_map[Counters::NGenLepInEtaPtAcc].first  = (float)n_gen_lep_in_eta_pt_acc;
   variables_map[Counters::NGenLepInEtaPtAcc].second = 1;

   variables_map[Counters::NGenLepNotInEtaPtAcc].first  = (float)(n_gen_lep - n_gen_lep_in_eta_pt_acc);
   variables_map[Counters::NGenLepNotInEtaPtAcc].second = 1;

   variables_map[Counters::NGenHLepNotInEtaPtAcc].first  = (float)(n_gen_H_lep - n_gen_H_lep_in_eta_pt_acc);
   variables_map[Counters::NGenHLepNotInEtaPtAcc].second = 1;

   variables_map[Counters::NGenAssocLepNotInEtaPtAcc].first  = (float)(n_gen_assoc_lep - n_gen_assoc_lep_in_eta_pt_acc);
   variables_map[Counters::NGenAssocLepNotInEtaPtAcc].second = 1;

   variables_map[Counters::NGenLepMinusNGoodLep].first  = (float)(n_gen_lep - (4 + nExtraLep));
   variables_map[Counters::NGenLepMinusNGoodLep].second = 1;

   variables_map[Counters::NGenLepInEtaPtAccMinusNGoodLep].first  = (float)(n_gen_lep_in_eta_pt_acc - (4 + nExtraLep));
   variables_map[Counters::NGenLepInEtaPtAccMinusNGoodLep].second = 1;

   variables_map[Counters::NExtraLep].first  = (float)nExtraLep;
   variables_map[Counters::NExtraLep].second = 1;

   variables_map[Counters::NExtraZ].first  = (float)nExtraZ;
   variables_map[Counters::NExtraZ].second = 1;

   variables_map[Counters::NJets].first  = (float)nCleanedJets;
   variables_map[Counters::NJets].second = 1;

   variables_map[Counters::NBtaggedJets].first  = (float)nCleanedJetsPt30BTagged;
   variables_map[Counters::NBtaggedJets].second = 1;

   variables_map[Counters::MET].first  = PFMET;
   variables_map[Counters::MET].second = 1;


//   cout << Counters::Pvbf1j << " " << variables_map[Counters::Pvbf1j].first << " " << variables_map[Counters::Pvbf1j].second << endl;

   for ( map<Counters::variable, pair<float, bool>>::iterator it = variables_map.begin(); it != variables_map.end(); it++ )
   {
      cout << it->second.second << endl;
      if ( !it->second.second ) continue;
      histograms->FillVariables( it->second.first, event_weight_, it->first, current_process_, reco_ch_1, reco_ch_2 );

//      cout << it->second.first << " => " << it->second.second << endl;
   }





//      for(int v=0; v<nVariables; v++){
//   if(!varPassCut[v]) continue;
//   hBCInSR[v][currentProcess][rc]->Fill(varVal[v],eventWeight);
//   hBCInSR[v][currentProcess][rc2]->Fill(varVal[v],eventWeight);
//      }
//
//      Bool_t varPairPassCut[nVariables] = {
//   1,1,nJets>=2,
//      };
//      for(int v2=0; v2<nVarPairs; v2++){
//   if(!varPairPassCut[v2]) continue;
//   h2DBCInSR[v2][currentProcess][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//   h2DBCInSR[v2][currentProcess][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//   if(!varPairRef[v2][2]) continue;
//   h2DBCInSRDecays[v2][currentProcess][0][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//   h2DBCInSRDecays[v2][currentProcess][0][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//   if(nGenHLep==4){
//     h2DBCInSRDecays[v2][currentProcess][1][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//     h2DBCInSRDecays[v2][currentProcess][1][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//     if(currentMatchHLepsStatus==0){
//       h2DBCInSRDecays[v2][currentProcess][3][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//       h2DBCInSRDecays[v2][currentProcess][3][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//     }else if(currentMatchHLepsStatus<5){
//       h2DBCInSRDecays[v2][currentProcess][4][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//       h2DBCInSRDecays[v2][currentProcess][4][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//     }
//   }else{
//     h2DBCInSRDecays[v2][currentProcess][2][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//     h2DBCInSRDecays[v2][currentProcess][2][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
//   }
//      }
//
//      nbWithBCInSRMatchHLeps[currentMatchHLepsStatus][currentProcess][rc]++;
//      nbWithBCInSRMatchHLeps[currentMatchHLepsStatus][currentProcess][rc2]++;
//      for(int v=0; v<nVariables; v++){
//   hBCInSRMatchHLeps[v][currentMatchHLepsStatus][currentProcess][rc]->Fill(varVal[v],eventWeight);
//   hBCInSRMatchHLeps[v][currentMatchHLepsStatus][currentProcess][rc2]->Fill(varVal[v],eventWeight);
//      }
//
//      nbWithBCInSRMatchAllLeps[currentMatchAllLepsStatus][currentProcess][rc]++;
//      nbWithBCInSRMatchAllLeps[currentMatchAllLepsStatus][currentProcess][rc2]++;
//      for(int v=0; v<nVariables; v++){
//   hBCInSRMatchAllLeps[v][currentMatchAllLepsStatus][currentProcess][rc]->Fill(varVal[v],eventWeight);
//   hBCInSRMatchAllLeps[v][currentMatchAllLepsStatus][currentProcess][rc2]->Fill(varVal[v],eventWeight);
//      }
//   
//      if(currentProcess==WH){
//   nbWithBCInSRMatchWH[currentMatchWHStatus][rc]++;
//   nbWithBCInSRMatchWH[currentMatchWHStatus][rc2]++;
//   for(int v=0; v<nVariables; v++){
//     hBCInSRMatchWH[v][currentMatchWHStatus][rc]->Fill(varVal[v],eventWeight);
//     hBCInSRMatchWH[v][currentMatchWHStatus][rc2]->Fill(varVal[v],eventWeight);
//   }
//      }
//      if(currentProcess==ZH){
//   nbWithBCInSRMatchZH[currentMatchZHStatus][rc]++;
//   nbWithBCInSRMatchZH[currentMatchZHStatus][rc2]++;
//   for(int v=0; v<nVariables; v++){
//     hBCInSRMatchZH[v][currentMatchZHStatus][rc]->Fill(varVal[v],eventWeight);
//     hBCInSRMatchZH[v][currentMatchZHStatus][rc2]->Fill(varVal[v],eventWeight);
//   }
//      }
//      if(currentProcess==ttH){
//   nbWithBCInSRMatchttH[currentMatchttHStatus][rc]++;
//   nbWithBCInSRMatchttH[currentMatchttHStatus][rc2]++;
//   for(int v=0; v<nVariables; v++){
//     hBCInSRMatchttH[v][currentMatchttHStatus][rc]->Fill(varVal[v],eventWeight);
//     hBCInSRMatchttH[v][currentMatchttHStatus][rc2]->Fill(varVal[v],eventWeight);
//   }
//      }
//
//      if(currentMatchAllLepsStatus==0){
//   nbWithBCInSRAll4LepRight[currentProcess][rc]++;
//   nbWithBCInSRAll4LepRight[currentProcess][rc2]++;
//      }
//
//      if(currentProcess==WH){
//   Int_t currentWDecay = -1;
//   if(nGenHLep==4 && nGenAssocLep==0) currentWDecay = 0;
//   else if(nGenHLep==4 && nGenAssocLep==1) currentWDecay = 1;
//   else cout<<"error"<<endl;
//   nbTotalWH[currentWDecay]++;
//   if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccWH[currentWDecay]++;
//   if(HLepsAreGood) nbHLepsAreGoodWH[currentWDecay]++;
//   if(currentMatchAllLepsStatus==0) nbAll4LepRightWH[currentWDecay]++;
//   nbZ1DaughtersFromHWH[currentWDecay][currentZ1MatchStatus]++;
//   nbZ2DaughtersFromHWH[currentWDecay][currentZ2MatchStatus]++;
//      }
//      if(currentProcess==ZH){
//   Int_t currentZDecay = -1;
//   if(nGenHLep==4 && nGenAssocLep==0) currentZDecay = 0;
//   else if(nGenHLep==4 && nGenAssocLep==2) currentZDecay = 1;
//   else if(nGenHLep==2 && nGenAssocLep==2) currentZDecay = 2;
//   else cout<<"error"<<endl;
//   nbTotalZH[currentZDecay]++;
//   if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccZH[currentZDecay]++;
//   if(HLepsAreGood) nbHLepsAreGoodZH[currentZDecay]++;
//   if(currentMatchAllLepsStatus==0) nbAll4LepRightZH[currentZDecay]++;
//   nbZ1DaughtersFromHZH[currentZDecay][currentZ1MatchStatus]++;
//   nbZ2DaughtersFromHZH[currentZDecay][currentZ2MatchStatus]++;
//      }
//      if(currentProcess==ttH){
//   Int_t currentttDecay = -1;
//   if(nGenHLep==4 && nGenAssocLep==0) currentttDecay = 0;
//   else if(nGenHLep==4 && nGenAssocLep==1) currentttDecay = 1;
//   else if(nGenHLep==4 && nGenAssocLep==2) currentttDecay = 2;
//   else if(nGenHLep==2 && nGenAssocLep==2) currentttDecay = 3;
//   else cout<<"error"<<endl;
//   nbTotalttH[currentttDecay]++;
//   if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccttH[currentttDecay]++;
//   if(HLepsAreGood) nbHLepsAreGoodttH[currentttDecay]++;
//   if(currentMatchAllLepsStatus==0) nbAll4LepRightttH[currentttDecay]++;
//   nbZ1DaughtersFromHttH[currentttDecay][currentZ1MatchStatus]++;
//   nbZ2DaughtersFromHttH[currentttDecay][currentZ2MatchStatus]++;
//      }
//
//      Bool_t rocPassCut[nROCs] = {
//   nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,
//      };
//      for(int ro=0; ro<nROCs; ro++){
//   if(!rocPassCut[ro]) continue;
//   if(currentProcess==rocRef[ro][1])
//     hRocSgnl[ro]->Fill(varVal[rocRef[ro][0]],eventWeight);
//   else if(currentProcess==rocRef[ro][2])
//     hRocBkgd[ro]->Fill(varVal[rocRef[ro][0]],eventWeight);
//      }
//      Bool_t evwPassCut[nEVWs] = {
//   nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nExtraLep>=1,nExtraZ>=1,
//      };
//      for(int e=0; e<nEVWs; e++){
//   if(!evwPassCut[e]) continue;
//     hEvw[e][currentProcess]->Fill(varVal[evwRef[e][0]],eventWeight);
//      }



}
//====================================



//==========================================
void Categorisation::ResetPerEventStuff()
{
   // Clean per-event counters
   n_gen_H_lep               = 0;
   n_gen_H_lep_in_eta_acc    = 0;
   n_gen_H_lep_in_pt_acc     = 0;
   n_gen_H_lep_in_eta_pt_acc = 0;
   n_gen_H_ele = 0;
   n_gen_H_mu  = 0;
   n_gen_H_tau = 0;
   n_gen_H_LEP = 0; // including tau

   n_gen_assoc_lep = 0;
   n_gen_assoc_lep_in_eta_acc    = 0;
   n_gen_assoc_lep_in_pt_acc     = 0;
   n_gen_assoc_lep_in_eta_pt_acc = 0;
   n_gen_assoc_ele = 0;
   n_gen_assoc_mu  = 0;
   n_gen_assoc_tau = 0;
   n_gen_assoc_LEP = 0; // including tau
   
   n_gen_LEP_plus  = 0; // including tau
   n_gen_LEP_minus = 0; // including tau
   
   n_gen_lep = 0;
   n_gen_lep_in_eta_acc    = 0;
   n_gen_lep_in_pt_acc     = 0;
   n_gen_lep_in_eta_pt_acc = 0;
   n_gen_ele = 0;
   n_gen_mu  = 0;
   n_gen_tau = 0;
   n_gen_LEP = 0; // including tau
   
   // Clean per-event maps
   counter_map.clear();
   
   // Clean per-event vectors
   gen_H_lep_id_.clear();
   gen_H_lep_pt_.clear();
   gen_H_lep_eta_.clear();
   gen_H_lep_phi_.clear();
   gen_assoc_lep_id_.clear();
   gen_assoc_lep_pt_.clear();
   gen_assoc_lep_eta_.clear();
   gen_assoc_lep_phi_.clear();
   
   sorted_gen_H_lep_pt_.clear();
   sorted_gen_H_lep_abs_eta_.clear();
}
//==========================================



//=======================================================
void Categorisation::SaveHistograms( TString file_name )
{
   histograms->SaveHistograms(file_name);
}
//=======================================================
