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
   
   m4l_min_ = 118.;
   m4l_max_ = 130.;
   
   gen_ch_1 = -999;
   
   reco_ch_1 = Counters::reco_ch_1_def;
   reco_ch_2 = Counters::reco_ch_2_def;
   reco_ch_3 = Counters::reco_ch_3_def;
   
   histograms = new Histograms(lumi_);
   
   signal_region = ZZsel >= 90;
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
   
   if (fChain == 0) return;

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
   bool found_matching_ambiguity = false;

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

   int n_ones = 0;
   int n_ones_H_lep = 0;
   int n_ones_assoc_lep = 0;
   
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
