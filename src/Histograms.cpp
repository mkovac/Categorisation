#include "Histograms.h"

using namespace std;

//Constructor
//==================================
Histograms::Histograms( float lumi)
{
   lumi_ = lumi;
   
   s_process_.push_back("Data");
   s_process_.push_back("H125");
   s_process_.push_back("H125ggH");
   s_process_.push_back("H125VBF");
   s_process_.push_back("H125WH");
   s_process_.push_back("H125ZH");
   s_process_.push_back("H125bbH");
   s_process_.push_back("H125ttH");
   s_process_.push_back("qqZZ");
   s_process_.push_back("ggZZ");
   s_process_.push_back("ttbar");
   
   
   // Sort
   s_sort_.push_back("first");
   s_sort_.push_back("second");
   s_sort_.push_back("third");
   s_sort_.push_back("fourth");
   
   
   // Generated channels
   s_gen_ch_.push_back("gen_ch_4mu");
   s_gen_ch_.push_back("gen_ch_4e");
   s_gen_ch_.push_back("gen_ch_2e2mu");
   s_gen_ch_.push_back("gen_ch_4tau");
   s_gen_ch_.push_back("gen_ch_2e2tau");
   s_gen_ch_.push_back("gen_ch_2mu2tau");
   s_gen_ch_.push_back("gen_ch_other");
   s_gen_ch_.push_back("gen_ch_e_mu");
   s_gen_ch_.push_back("gen_ch_tau");
   s_gen_ch_.push_back("gen_ch_2");
   s_gen_ch_.push_back("gen_ch_3");
   
   
   // Reconstructed channels
   s_reco_ch_.push_back("reco_ch_4mu");
   s_reco_ch_.push_back("reco_ch_4e");
   s_reco_ch_.push_back("reco_ch_2e2mu");
   s_reco_ch_.push_back("reco_ch_all");
   s_reco_ch_.push_back("reco_ch_1_def");
   s_reco_ch_.push_back("reco_ch_2_def");
   s_reco_ch_.push_back("reco_ch_3_def");
   
   // Decays
   s_variable_decay_.push_back("H_to_any");
   s_variable_decay_.push_back("H_to_4l");
   s_variable_decay_.push_back("H_to_2l2X");
   s_variable_decay_.push_back("H_to_4l_4_from_H");
   s_variable_decay_.push_back("H_to_4l_not_4_from_H");

   // Match H leptons status
   s_H_lep_match_status_.push_back("H_4");
   s_H_lep_match_status_.push_back("H_3");
   s_H_lep_match_status_.push_back("H_2");
   s_H_lep_match_status_.push_back("H_1");
   s_H_lep_match_status_.push_back("H_0");
   s_H_lep_match_status_.push_back("ambiguity");

   // Match all leptons status
   s_all_lep_match_status_.push_back("all_40");
   s_all_lep_match_status_.push_back("all_31");
   s_all_lep_match_status_.push_back("all_22");
   s_all_lep_match_status_.push_back("all_le_4");
   s_all_lep_match_status_.push_back("ambiguity");
   
   // Match WH leptons status
   s_WH_lep_match_status_.push_back("WH_40_40");
   s_WH_lep_match_status_.push_back("WH_40_41");
   s_WH_lep_match_status_.push_back("WH_31_41");
   s_WH_lep_match_status_.push_back("WH_le_4");
   s_WH_lep_match_status_.push_back("ambiguity");
   
   // Match ZH leptons status
   s_ZH_lep_match_status_.push_back("ZH_40_40");
   s_ZH_lep_match_status_.push_back("ZH_40_42");
   s_ZH_lep_match_status_.push_back("ZH_31_42");
   s_ZH_lep_match_status_.push_back("ZH_22_42");
   s_ZH_lep_match_status_.push_back("ZH_22_22");
   s_ZH_lep_match_status_.push_back("ZH_le_4");
   s_ZH_lep_match_status_.push_back("ambiguity");
   
   // Match ttH leptons status
   s_ttH_lep_match_status_.push_back("ttH_40_40");
   s_ttH_lep_match_status_.push_back("ttH_40_41");
   s_ttH_lep_match_status_.push_back("ttH_31_41");
   s_ttH_lep_match_status_.push_back("ttH_40_42");
   s_ttH_lep_match_status_.push_back("ttH_31_42");
   s_ttH_lep_match_status_.push_back("ttH_22_42");
   s_ttH_lep_match_status_.push_back("ttH_22_22");
   s_ttH_lep_match_status_.push_back("ttH_le_4");
   s_ttH_lep_match_status_.push_back("ambiguity");
   
   //Prepare variable pairs
   var_pair.PrepareVarPair("M4l_vs_D_kin_bkg","m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", 100, 50, 850, 20, 0, 1);
   var_pair.PrepareVarPair("MZ2_vs_D_kin_bkg", "m_{Z_{2}} (GeV)", "D_{bkg}^{kin}", 75, 0, 150, 20, 0, 1);
   var_pair.PrepareVarPair("D2jVbfHjj_vs_D2jqg", "D_{VBF-2j}^{ME}", "D_{2jets}^{q/g}", 51, 0, 1.02, 25, -0.04, 0.96);

   //Prepare variables
   variable.PrepareVar("M4l", "m_{4#font[12]{l}} (GeV)", 100, 50, 85);
   variable.PrepareVar("M4l2", "m_{4#font[12]{l}} (GeV)", 70, 105, 140);
   variable.PrepareVar("MZ1", "m_{Z_{1}} (GeV)", 75, 0, 150);
   variable.PrepareVar("MZ2", "m_{Z_{2}} (GeV)", 75, 0, 150);
   variable.PrepareVar("Dkinbkg", "D_{bkg}^{kin}", 20, 0, 1);
   variable.PrepareVar("DjetFisher", "D^{Fisher}", 50, 0, 2);
   variable.PrepareVar("Pvbf", "P_{VBF}^{MELA}", 50, -2, 98);
   variable.PrepareVar("Phjj", "P_{ggH+2j}^{MELA}", 50, -5, 245);
   variable.PrepareVar("Pvbf1j", "P_{VBF 1j}^{MELA}", 50, -10, 490);
   variable.PrepareVar("Phj", "P_{ggH+1j}^{MELA}", 50, -10, 490);
   variable.PrepareVar("Pwhhadr", "P_{WH-h}^{MELA}", 50, -10, 490);
   variable.PrepareVar("Pzhhadr", "P_{ZH-h}^{MELA}", 50, -10, 490);
   variable.PrepareVar("Pwhlept", "P_{WH-l}^{MELA}", 50, -1, 49);
   variable.PrepareVar("Pzhlept", "P_{ZH-l}^{MELA}", 50, -1, 49);
   variable.PrepareVar("D2jVbfHjj", "D_{VBF-2j}^{ME}", 51, 0, 1.02);
   variable.PrepareVar("D1jVbfHj", "D_{VBF-1j}^{ME}", 51, 0, 1.02);
   variable.PrepareVar("D2jWHHadrHjj", "D_{WH-hadr.}^{ME}", 51, 0, 1.02);
   variable.PrepareVar("D2jZHHadrHjj", "D_{ZH-hadr.}^{ME}", 51, 0, 1.02);
   variable.PrepareVar("Pqj1", "P_{q}(j_{1})", 25, -0.04, 0.96);
   variable.PrepareVar("Pgj1", "P_{g}(j_{1})", 25, -0.04, 0.96);
   variable.PrepareVar("Pqj1Pqj2", "P_{q}(j_{1})*P_{q}(j_{2})", 25, -0.04, 0.96);
   variable.PrepareVar("Pgj1Pgj2", "P_{g}(j_{1})*P_{g}(j_{2})", 25, -0.04, 0.96);
   variable.PrepareVar("D2jqg", "D_{2jets}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("Dqgj1Dqgj2", "D_{1jet}^{q/g}(j_{1})*D_{1jet}^{q/g}(j_{2})", 51, 0, 1.02);
   variable.PrepareVar("Pqj1VbfTopo", "P_{q}(j_{1})", 25, -0.04, 0.96);
   variable.PrepareVar("Pgj1VbfTopo", "P_{g}(j_{1})", 25, -0.04, 0.96);
   variable.PrepareVar("Pqj1Pqj2VbfTopo", "P_{q}(j_{1})*P_{q}(j_{2})", 25, -0.04, 0.96);
   variable.PrepareVar("Pgj1Pgj2VbfTopo", "P_{g}(j_{1})*P_{g}(j_{2})", 25, -0.04, 0.96);
   variable.PrepareVar("D2jqgVbfTopo", "D_{2jets}^{q/g}", 26, 0, 1.04);
   variable.PrepareVar("Pq", "P_{q}", 25, -0.04, 0.96);
   variable.PrepareVar("Pg", "P_{g}", 25, -0.04, 0.96);
   variable.PrepareVar("D1jqg", "D_{1jet}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQGVbfHjj", "D_{VBF-2j}^{ME q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaD2jQGVbfHjj", "D_{VBF-2j}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaQGVbfHj", "D_{VBF-1j}^{ME q/g}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaD1jQGVbfHj", "D_{VBF-1j}^{ME}*D_{1jet}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQGWHHadrHjj", "D_{WH-hadr.}^{ME q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaD2jQGWHHadrHjj", "D_{WH-hadr.}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQGZHHadrHjj", "D_{ZH-hadr.}^{ME q/g}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaD2jQGZHHadrHjj", "D_{ZH-hadr.}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02);
   variable.PrepareVar("RatioPvbfPhjj", "P_{VBF}^{MELA} / P_{ggH+2j}^{MELA}", 50, 0, 5);
   variable.PrepareVar("RatioPqj1Pqj2Pgj1Pgj2", "[P_{q}(j_{1})*P_{q}(j_{2})] / [P_{g}(j_{1})*P_{g}(j_{2})]", 50, 0, 5);
   variable.PrepareVar("RatioPvbfPqj1Pqj2PhjjPgj1Pgj2", "[P_{VBF}^{MELA}*P_{q}(j_{1})*P_{q}(j_{2})] / [P_{ggH+2j}^{MELA}*P_{g}(j_{1})*P_{g}(j_{2})]", 50, 0, 5);
   variable.PrepareVar("D2jMelaExpQGVbfHjj", "D_{VBF-2j}^{ME exp(q/g)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaSqQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^2}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaSqrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/2)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaCbrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/3)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQrrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/4)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQnrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/5)}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaSqrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/2)}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaCbrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/3)}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaQrrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/4)}", 51, 0, 1.02);
   variable.PrepareVar("D1jMelaQnrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/5)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaSqrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/2)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaCbrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/3)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQrrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/4)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQnrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/5)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaSqrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/2)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaCbrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/3)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQrrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/4)}", 51, 0, 1.02);
   variable.PrepareVar("D2jMelaQnrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/5)}", 51, 0, 1.02);
   variable.PrepareVar("Pt4l", "p_{T}^{4l} (GeV)", 50, 0, 500);
   variable.PrepareVar("NGenLep", "# gen leptons", 7, 0, 7);
   variable.PrepareVar("NGenLepInEtaPtAcc", "# gen leptons in acceptance", 7, 0, 7);
   variable.PrepareVar("NGenLepNotInEtaPtAcc", "# gen leptons not in acceptance", 7, 0, 7);
   variable.PrepareVar("NGenHLepNotInEtaPtAcc", "# gen leptons from H not in acceptance", 5, 0, 5);
   variable.PrepareVar("NGenAssocLepNotInEtaPtAcc", "# gen associated leptons not in acceptance", 3, 0, 3);
   variable.PrepareVar("NGenLepMinusNGoodLep", "# gen leptons - # good leptons", 6, -3, 3);
   variable.PrepareVar("NGenLepInEtaPtAccMinusNGoodLep", "# gen leptons in acceptance - # good leptons", 6, -3, 3);
   variable.PrepareVar("NExtraLep", "number of additional leptons", 6, 0, 6);
   variable.PrepareVar("NExtraZ", "number of additional #font[12]{l}^{+}#font[12]{l}^{-} pairs", 4, 0, 4);
   variable.PrepareVar("NJets", "number of selected jets", 15, 0, 15);
   variable.PrepareVar("NBtaggedJets", "number of selected b-tagged jets", 6, 0, 6);
   variable.PrepareVar("MET", "E_{T}^{miss.}", 100, 0, 200);

   
   // Prepare histograms
   for ( int i_proc = 0; i_proc < Counters::num_of_processes; i_proc++ )
   {
      histo_name_ = "h_gen_H_pt_" + s_process_.at(i_proc);
      histo_label_ = ";" + Variables::pt().var_X_label + ";" + Variables::pt().var_Y_label;
      h_gen_H_pt_[i_proc] = new TH1F(histo_name_, histo_label_, 100, 0, 100);
      
      histo_name_ = "h_gen_H_eta_" + s_process_.at(i_proc);
      histo_label_ = ";" + Variables::eta().var_X_label + ";" + Variables::eta().var_Y_label;
      h_gen_H_eta_[i_proc] = new TH1F(histo_name_, histo_label_, 100, 0, 5);
      
      histo_name_ = "h_gen_H_eta_vs_pt_" + s_process_.at(i_proc);
      histo_label_ = ";" + Variables::pt().var_X_label + ";" + Variables::eta().var_X_label;
      h_gen_H_eta_vs_pt_[i_proc] = new TH2F(histo_name_, histo_label_, 20, 0, 100, 20, 0, 5);
      
      for ( int i_rc = 0; i_rc < Counters::num_of_reco_ch; i_rc++ )
      {
         for ( int i_var = 0; i_var < Counters::num_of_vars; i_var++ )
         {
            histo_label_ = ";" + variable.vec_var[i_var].label + ";" + "# of events";

            histo_name_ = "h_bc_in_sig_reg_" + variable.vec_var[i_var].name + "_" +  s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
            h_bc_in_sig_reg_[i_var][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                        variable.vec_var[i_var].x_min,
                                                                                        variable.vec_var[i_var].x_max);
            
            for ( int i_ms = 0; i_ms < Counters::num_of_H_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_H_leps_" + variable.vec_var[i_var].name + "_" + s_H_lep_match_status_.at(i_ms) + "_"
                                                             + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_H_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                                              variable.vec_var[i_var].x_min,
                                                                                                              variable.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_all_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_all_leps_" + variable.vec_var[i_var].name + "_" + s_all_lep_match_status_.at(i_ms) + "_"
                                                               + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_all_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                                                variable.vec_var[i_var].x_min,
                                                                                                                variable.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_WH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_WH_leps_" + variable.vec_var[i_var].name + "_" + s_WH_lep_match_status_.at(i_ms) + "_"
                                                              + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_WH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                                               variable.vec_var[i_var].x_min,
                                                                                                               variable.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_ZH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_ZH_leps_" + variable.vec_var[i_var].name + "_" + s_ZH_lep_match_status_.at(i_ms) + "_"
                                                              + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_ZH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                                               variable.vec_var[i_var].x_min,
                                                                                                               variable.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_ttH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_ttH_leps_" + variable.vec_var[i_var].name + "_" + s_ttH_lep_match_status_.at(i_ms) + "_"
                                                               + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_ttH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, variable.vec_var[i_var].n_bins,
                                                                                                                variable.vec_var[i_var].x_min,
                                                                                                                variable.vec_var[i_var].x_max);
            }
         } // end i_var
         
         // Variable pairs histograms
         for ( int i_var_pair = 0; i_var_pair < Counters::num_of_variable_pairs; i_var_pair++ )
         {
            histo_name_ = "h_bc_in_sig_reg_2D_" + var_pair.vec_var_pair[i_var_pair].name + "_" +  s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
            histo_label_ = ";" + var_pair.vec_var_pair[i_var_pair].x_label + ";" + var_pair.vec_var_pair[i_var_pair].y_label;
//            if ( histo_name_.Contains("ggH") ) cout << histo_name_ << endl;
            h_bc_in_sig_reg_2D_[i_var_pair][i_proc][i_rc] = new TH2F(histo_name_, histo_label_, var_pair.vec_var_pair[i_var_pair].x_bins,
                                                                                                var_pair.vec_var_pair[i_var_pair].x_min,
                                                                                                var_pair.vec_var_pair[i_var_pair].x_max,
                                                                                                var_pair.vec_var_pair[i_var_pair].y_bins,
                                                                                                var_pair.vec_var_pair[i_var_pair].y_min,
                                                                                                var_pair.vec_var_pair[i_var_pair].y_max);
            
            for ( int i_decay = 0; i_decay < 5; i_decay++ )
            {
               histo_name_ = "h_bc_in_sig_reg_2D_decays_" + var_pair.vec_var_pair[i_var_pair].name + "_" +  s_process_.at(i_proc) + "_"
                                                          + s_variable_decay_.at(i_decay) + "_" + s_reco_ch_.at(i_rc);
               histo_label_ = ";" + var_pair.vec_var_pair[i_var_pair].x_label + ";" + var_pair.vec_var_pair[i_var_pair].y_label;
//               if ( histo_name_.Contains("ggH") ) cout << histo_name_ << endl;
               h_bc_in_sig_reg_2D_decays_[i_var_pair][i_proc][i_decay][i_rc] = new TH2F(histo_name_, histo_label_, var_pair.vec_var_pair[i_var_pair].x_bins,
                                                                                                                   var_pair.vec_var_pair[i_var_pair].x_min,
                                                                                                                   var_pair.vec_var_pair[i_var_pair].x_max,
                                                                                                                   var_pair.vec_var_pair[i_var_pair].y_bins,
                                                                                                                   var_pair.vec_var_pair[i_var_pair].y_min,
                                                                                                                   var_pair.vec_var_pair[i_var_pair].y_max);
            } // end i_decay
         } // end i_var_pair
      } // end i_rc
      
   
      for ( int i_gc = 0; i_gc < Counters::num_of_gen_ch; i_gc++ )
      {
         histo_name_ = "h_num_reco_H_ele_in_eta_pt_acc_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc);
         histo_label_ = ";" + Variables::n_ele().var_X_label + ";" + Variables::n_ele().var_Y_label;
         h_num_reco_H_ele_in_eta_pt_acc_[i_proc][i_gc] = new TH1F(histo_name_, histo_label_, Variables::n_ele().var_N_bin,
                                                                  Variables::n_ele().var_min, Variables::n_ele().var_max);
         
         histo_name_ = "h_num_reco_H_mu_in_eta_pt_acc_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc);
         histo_label_ = ";" + Variables::n_mu().var_X_label + ";" + Variables::n_mu().var_Y_label;
         h_num_reco_H_mu_in_eta_pt_acc_[i_proc][i_gc] = new TH1F(histo_name_, histo_label_, Variables::n_mu().var_N_bin,
                                                                  Variables::n_mu().var_min, Variables::n_mu().var_max);

         histo_name_ = "h_num_reco_H_lep_in_eta_pt_acc_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc);
         histo_label_ = ";" + Variables::n_lep().var_X_label + ";" + Variables::n_lep().var_Y_label;
         h_num_reco_H_lep_in_eta_pt_acc_[i_proc][i_gc] = new TH1F(histo_name_, histo_label_, Variables::n_lep().var_N_bin,
                                                                  Variables::n_lep().var_min, Variables::n_lep().var_max);

         for ( int i_sort = 0; i_sort < Counters::num_of_sorted_objects; i_sort++ )
         {
            histo_name_ = "h_pt_gen_H_lep_in_eta_acc_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc) + "_" + s_sort_.at(i_sort);
            histo_label_ = ";" + Variables::pt().var_X_label + ";" + Variables::pt().var_Y_label;
            h_pt_gen_H_lep_in_eta_acc_[i_proc][i_gc][i_sort] = new TH1F(histo_name_, histo_label_, Variables::pt().var_N_bin,
                                                                        Variables::pt().var_min, Variables::pt().var_max);
            
            histo_name_ = "h_pt_reco_bc_in_sig_reg_and_pass_triger_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc) + "_" + s_sort_.at(i_sort);
            histo_label_ = ";" + Variables::pt().var_X_label + ";" + Variables::pt().var_Y_label;
            h_pt_reco_bc_in_sig_reg_and_pass_triger_[i_proc][i_gc][i_sort] = new TH1F(histo_name_, histo_label_, Variables::pt().var_N_bin,
                                                                        Variables::pt().var_min, Variables::pt().var_max);
            
            
            histo_name_ = "h_eta_gen_H_lep_in_pt_acc_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc) + "_" + s_sort_.at(i_sort);
            histo_label_ = ";" + Variables::eta().var_X_label + ";" + Variables::eta().var_Y_label;
            h_eta_gen_H_lep_in_pt_acc_[i_proc][i_gc][i_sort] = new TH1F(histo_name_, histo_label_, Variables::eta().var_N_bin,
                                                                        Variables::eta().var_min, Variables::eta().var_max);
            
            histo_name_ = "h_eta_reco_bc_in_sig_reg_and_pass_triger_" + s_process_.at(i_proc) + "_" + s_gen_ch_.at(i_gc) + "_" + s_sort_.at(i_sort);
            histo_label_ = ";" + Variables::eta().var_X_label + ";" + Variables::eta().var_Y_label;
            h_eta_reco_bc_in_sig_reg_and_pass_triger_[i_proc][i_gc][i_sort] = new TH1F(histo_name_, histo_label_, Variables::eta().var_N_bin,
                                                                        Variables::eta().var_min, Variables::eta().var_max);
         } // end i_sort
      } // end i_gc
   } // end i_proc
   
   
//   for ( int i_roc = 0; i_roc < Counters::num_of_ROCs; i_roc++ )
//   {
//      histo_name_ = "h_ROC_sig_" + s_ROC_.at(i_roc);
//      histo_label_ = ";" + Variables::eta().var_X_label + ";" + Variables::eta().var_Y_label;
//      h_ROC_sig_[i_roc] = new TH1F(histo_name_, histo_label_, Variables::eta().var_N_bin,
//                                                                        Variables::eta().var_min, Variables::eta().var_max);
//      
//   }
   
   
   
   
   
}
//==================================



//=======================
Histograms::~Histograms()
{
}
//=======================



//=================================================================================================
void Histograms::FillPt( float pT, float weight, int proc, int gc_1, int gc_2, int gc_3, int sort )
{
   h_pt_gen_H_lep_in_eta_acc_[proc][gc_1][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
   h_pt_gen_H_lep_in_eta_acc_[proc][gc_2][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
   h_pt_gen_H_lep_in_eta_acc_[proc][gc_3][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillPtReco( float pT, float weight, int proc, int gc_1, int gc_2, int gc_3, int sort )
{
   h_pt_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_1][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
   h_pt_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_2][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
   h_pt_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_3][sort]->Fill(pT, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillAbsEta( float eta, float weight, int proc, int gc_1, int gc_2, int gc_3, int sort )
{
   h_eta_gen_H_lep_in_pt_acc_[proc][gc_1][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
   h_eta_gen_H_lep_in_pt_acc_[proc][gc_2][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
   h_eta_gen_H_lep_in_pt_acc_[proc][gc_3][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillAbsEtaReco( float eta, float weight, int proc, int gc_1, int gc_2, int gc_3, int sort )
{
   h_eta_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_1][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
   h_eta_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_2][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
   h_eta_reco_bc_in_sig_reg_and_pass_triger_[proc][gc_3][sort]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillRecoEleN( int num, float weight, int proc, int gc_1, int gc_2, int gc_3 )
{
   h_num_reco_H_ele_in_eta_pt_acc_[proc][gc_1]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_ele_in_eta_pt_acc_[proc][gc_2]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_ele_in_eta_pt_acc_[proc][gc_3]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillRecoMuN( int num, float weight, int proc, int gc_1, int gc_2, int gc_3 )
{
   h_num_reco_H_mu_in_eta_pt_acc_[proc][gc_1]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_mu_in_eta_pt_acc_[proc][gc_2]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_mu_in_eta_pt_acc_[proc][gc_3]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillRecoLepN( int num, float weight, int proc, int gc_1, int gc_2, int gc_3 )
{
   h_num_reco_H_lep_in_eta_pt_acc_[proc][gc_1]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_lep_in_eta_pt_acc_[proc][gc_2]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
   h_num_reco_H_lep_in_eta_pt_acc_[proc][gc_3]->Fill(num, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillPtEtaH( float pt, float eta, float weight, int proc )
{
   h_gen_H_pt_[proc]->Fill(pt, (proc == Counters::AllData) ? 1. : weight);
   h_gen_H_eta_[proc]->Fill(eta, (proc == Counters::AllData) ? 1. : weight);
   h_gen_H_eta_vs_pt_[proc]->Fill(pt, eta, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillVariables( float var_value, float weight, int variable, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_[variable][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_[variable][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillMatchLepsH( float var_value, float weight, int variable, int match_status, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_match_H_leps_[variable][match_status][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_match_H_leps_[variable][match_status][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillMatchLepsAll( float var_value, float weight, int variable, int match_status, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_match_all_leps_[variable][match_status][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_match_all_leps_[variable][match_status][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillMatchLepsWH( float var_value, float weight, int variable, int match_status, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_match_WH_leps_[variable][match_status][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_match_WH_leps_[variable][match_status][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillMatchLepsZH( float var_value, float weight, int variable, int match_status, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_match_ZH_leps_[variable][match_status][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_match_ZH_leps_[variable][match_status][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillMatchLepsttH( float var_value, float weight, int variable, int match_status, int proc, int rc_1, int rc_2 )
{
   h_bc_in_sig_reg_match_ttH_leps_[variable][match_status][proc][rc_1]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_match_ttH_leps_[variable][match_status][proc][rc_2]->Fill(var_value, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillVariablePairs( float var_value_x, float var_value_y, float weight, int var_pair, int proc, int rc_1, int rc_2 )
{

//   cout << var_value_x << "   " << var_value_y << endl;
   h_bc_in_sig_reg_2D_[var_pair][proc][rc_1]->Fill(var_value_x, var_value_y, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_2D_[var_pair][proc][rc_2]->Fill(var_value_x, var_value_y, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//=================================================================================================
void Histograms::FillVariablePairsDecay( float var_value_x, float var_value_y, float weight, int var_pair, int proc, int decay, int rc_1, int rc_2 )
{

//   cout << var_value_x << "   " << var_value_y << endl;
   h_bc_in_sig_reg_2D_decays_[var_pair][proc][decay][rc_1]->Fill(var_value_x, var_value_y, (proc == Counters::AllData) ? 1. : weight);
   h_bc_in_sig_reg_2D_decays_[var_pair][proc][decay][rc_2]->Fill(var_value_x, var_value_y, (proc == Counters::AllData) ? 1. : weight);
}
//=================================================================================================



//==================================================
void Histograms::SaveHistograms( TString file_name )
{

   cout << "[INFO] Saving histograms to " << file_name << endl;
   
   TFile *f_out_histos = new TFile(file_name, "recreate");
   f_out_histos->cd();

   for ( int i_proc = 0; i_proc < Counters::num_of_processes; i_proc++ )
   {
      h_gen_H_pt_[i_proc]->Write();
      h_gen_H_eta_[i_proc]->Write();
      h_gen_H_eta_vs_pt_[i_proc]->Write();
   
     for ( int i_rc = 0; i_rc < Counters::num_of_reco_ch; i_rc++ )
      {
         for ( int i_var = 0; i_var < Counters::num_of_vars; i_var++ )
         {
            h_bc_in_sig_reg_[i_var][i_proc][i_rc]->Write();
            
            for ( int i_ms = 0; i_ms < Counters::num_of_H_lep_match_statuses; i_ms++ )
            {
               h_bc_in_sig_reg_match_H_leps_[i_var][i_ms][i_proc][i_rc]->Write();
            }
            for ( int i_ms = 0; i_ms < Counters::num_of_all_lep_match_statuses; i_ms++ )
            {
               h_bc_in_sig_reg_match_all_leps_[i_var][i_ms][i_proc][i_rc]->Write();
            }
            for ( int i_ms = 0; i_ms < Counters::num_of_WH_lep_match_statuses; i_ms++ )
            {
               h_bc_in_sig_reg_match_WH_leps_[i_var][i_ms][i_proc][i_rc]->Write();
            }
            for ( int i_ms = 0; i_ms < Counters::num_of_ZH_lep_match_statuses; i_ms++ )
            {
               h_bc_in_sig_reg_match_ZH_leps_[i_var][i_ms][i_proc][i_rc]->Write();
            }
            for ( int i_ms = 0; i_ms < Counters::num_of_ttH_lep_match_statuses; i_ms++ )
            {
               h_bc_in_sig_reg_match_ttH_leps_[i_var][i_ms][i_proc][i_rc]->Write();
            }
         } // end i_var
         
         for ( int i_var_pair = 0; i_var_pair < Counters::num_of_variable_pairs; i_var_pair++ )
         {
            h_bc_in_sig_reg_2D_[i_var_pair][i_proc][i_rc]->Write();
            
            for ( int i_decay = 0; i_decay < 5; i_decay++ )
            {
               h_bc_in_sig_reg_2D_decays_[i_var_pair][i_proc][i_decay][i_rc]->Write();
            } // end i_decay
         } // end i_var_pair
      } // end i_rc


      for ( int i_gc = 0; i_gc < Counters::num_of_gen_ch; i_gc++ )
      {
         h_num_reco_H_ele_in_eta_pt_acc_[i_proc][i_gc]->Write();
         h_num_reco_H_mu_in_eta_pt_acc_[i_proc][i_gc]->Write();
         h_num_reco_H_lep_in_eta_pt_acc_[i_proc][i_gc]->Write();

//         cout << Counters::num_of_sorted_objects << endl;

         for ( int i_sort = 0; i_sort < Counters::num_of_sorted_objects; i_sort++ )
         {
            h_pt_gen_H_lep_in_eta_acc_[i_proc][i_gc][i_sort]->Write();
            h_pt_reco_bc_in_sig_reg_and_pass_triger_[i_proc][i_gc][i_sort]->Write();
            h_eta_gen_H_lep_in_pt_acc_[i_proc][i_gc][i_sort]->Write();
            h_eta_reco_bc_in_sig_reg_and_pass_triger_[i_proc][i_gc][i_sort]->Write();
         }
      }
   }
   
   f_out_histos->Close();
   delete f_out_histos;
}
//==================================================



//=============================
//void Histograms::DeleteHistos()
//{
//   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
//   {
//      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
//      {
//         //=============
//         // M4l
//         //=============
//         delete histos_1D_ZX[Settings::M4lMain][i_fs][i_cat];
//         delete histos_1D_ZX_shape[Settings::M4lMain][i_fs][i_cat];
//         delete histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat];
//         delete histos_1D_ZX_shape[Settings::M4lMainZoomed][i_fs][i_cat];
//         delete histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat];
//         delete histos_1D_ZX_shape[Settings::M4lMainHighMass][i_fs][i_cat];
//
//
//      }
//   }
//
//}
//=============================




//=============================================
//void Histograms::GetHistos( TString file_name )
//{
//   TFile* histo_file = TFile::Open(file_name);
//
//   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
//   {
//      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
//      {
//         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
//         {
//            //=============
//            // M4l
//            //=============
//            _histo_name = "M4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
//            histos_1D[Settings::M4lMain][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
//
//            _histo_name = "M4l_zoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
//            histos_1D[Settings::M4lMainZoomed][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
//
//            _histo_name = "M4l_HighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
//            histos_1D[Settings::M4lMainHighMass][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
//
//      }
//   }
//}
//=============================================





//=========================================================================================================
//void Histograms::plot_1D_single( TString filename, TString variable_name, TString folder, int fs, int cat )
//{
//   int plot_index = SetPlotName( variable_name);
//
//   TCanvas *c;
//   if(variable_name == "M4lMain") c = new TCanvas(variable_name, variable_name, 650, 500);
//   else c = new TCanvas(variable_name, variable_name, 600, 600);
//
//   if ( GetVarLogX( variable_name) ) c->SetLogx();
//   if ( GetVarLogY( variable_name) ) c->SetLogy();
//
//
//   is_Djet_ = plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::D2jet_M4L118130 || plot_index == Settings::D1jet || plot_index == Settings::D2jet;
//   is_DVH_  = plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DWH || plot_index == Settings::DZH_M4L118130 || plot_index == Settings::DZH ||
//              plot_index == Settings::DVH_M4L118130 || plot_index == Settings::DVH;
//
//   if ( is_Djet_ )
//   {
//      histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetFillColor(Cosmetics::VBF().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125VH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//   }
//   else if ( is_DVH_ )
//   {
//      histos_1D[plot_index][fs][cat][Settings::H125VH]->SetFillColor(Cosmetics::VH().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//      histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
//   }
//   else
//   {
//      histos_1D[plot_index][fs][cat][Settings::H125]->SetFillColor(Cosmetics::Higgs_all().fill_color);
//   }
//
//   histos_1D[plot_index][fs][cat][Settings::qqZZ]->SetFillColor(Cosmetics::qqZZ().fill_color);
//   histos_1D[plot_index][fs][cat][Settings::ggZZ]->SetFillColor(Cosmetics::ggZZ().fill_color);
//   histos_1D[plot_index][fs][cat][Settings::qqZZ]->SetLineColor(Cosmetics::qqZZ().line_color);
//   histos_1D[plot_index][fs][cat][Settings::ggZZ]->SetLineColor(Cosmetics::ggZZ().line_color);
//
//   if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
//   {
//      histos_1D_ZX_shape[plot_index][fs][cat]->SetFillColor(Cosmetics::ZX().fill_color);
//      histos_1D_ZX_shape[plot_index][fs][cat]->SetLineColor(Cosmetics::ZX().line_color);
//   }
//   else
//   {
//      histos_1D_ZX[plot_index][fs][cat]->SetFillColor(Cosmetics::ZX().fill_color);
//      histos_1D_ZX[plot_index][fs][cat]->SetLineColor(Cosmetics::ZX().line_color);
//   }
//
//
//   histos_1D[plot_index][fs][cat][Settings::H125]->SetLineColor(Cosmetics::Higgs_all().line_color);
//   histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetLineColor(Cosmetics::Higgs_all().line_color);
//   histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetLineColor(Cosmetics::Higgs_all().line_color);
//   histos_1D[plot_index][fs][cat][Settings::H125VH]->SetLineColor(Cosmetics::Higgs_all().line_color);
//   histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetLineColor(Cosmetics::Higgs_all().line_color);
//
//   histos_1D[plot_index][fs][cat][Settings::AllData]->SetBinErrorOption(TH1::kPoisson);
//   histos_1D[plot_index][fs][cat][Settings::AllData]->SetLineColor(kBlack);
//
//
//   // THStack
//   THStack *stack = new THStack( "stack", "stack" );
//
//   if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
//   {
//      stack->Add(histos_1D_ZX_shape[plot_index][fs][cat]);
//   }
//   else
//   {
//      stack->Add(histos_1D_ZX[plot_index][fs][cat]);
//   }
//
//   stack->Add(histos_1D[plot_index][fs][cat][Settings::ggZZ]);
//   stack->Add(histos_1D[plot_index][fs][cat][Settings::qqZZ]);
//
//   if ( is_Djet_ )
//   {
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125VH]);
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125ttH]);
//      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125ggH]);
//      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125VBF]);
//   }
//   else if ( is_DVH_ )
//   {
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125VBF]);
//      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125ttH]);
//      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125ggH]);
//      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125VH]);
//   }
//   else
//   {
//      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125]);
//   }
//
//   stack->Draw("HIST");
//
//   float data_max = histos_1D[plot_index][fs][cat][Settings::AllData]->GetBinContent(histos_1D[plot_index][fs][cat][Settings::AllData]->GetMaximumBin());
//   float data_max_error = histos_1D[plot_index][fs][cat][Settings::AllData]->GetBinErrorUp(histos_1D[plot_index][fs][cat][Settings::AllData]->GetMaximumBin());
//
//   if ( GetVarLogY(variable_name) )
//   {
//      stack->SetMinimum(0.2);
//      stack->SetMaximum((data_max + data_max_error)*100);
//   }
//   else
//   {
//      stack->SetMinimum(1e-5);
//      stack->SetMaximum((data_max + data_max_error)*1.1);
//   }
//
//   if ( plot_index == Settings::M4lMain || plot_index == Settings::M4lMainZoomed )
//   {
//      if      ( fs == Settings::fs4e )    stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4e);
//      else if ( fs == Settings::fs4mu )   stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4mu);
//      else if ( fs == Settings::fs2e2mu ) stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_2e2mu);
//      else                                stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label);;
//   }
//   else
//   {
//      stack->GetXaxis()->SetTitle(histos_1D[plot_index][fs][cat][Settings::AllData]->GetXaxis()->GetTitle());
//   }
//
//   stack->GetYaxis()->SetTitle(histos_1D[plot_index][fs][cat][Settings::AllData]->GetYaxis()->GetTitle());
//
//   if ( (plot_index == Settings::M4lMainZoomed) || (plot_index == Settings::M4lMainHighMass) ) stack->GetXaxis()->SetNdivisions(1005);
//
//   histos_1D[plot_index][fs][cat][Settings::AllData]->Draw("SAME p E1 X0");
//
//
////=============
//// L E G E N D
////=============
//
//   TLegend *legend;
//
//   if(variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
//   {
//      legend  = CreateLegend("right", histos_1D[plot_index][fs][cat][Settings::AllData],
//                                      histos_1D[plot_index][fs][cat][Settings::H125],
//                                      histos_1D[plot_index][fs][cat][Settings::qqZZ],
//                                      histos_1D[plot_index][fs][cat][Settings::ggZZ],
//                                      histos_1D_ZX_shape[plot_index][fs][cat]);
//   }
//   else if ( plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::D1jet )
//   {
//      legend  = CreateLegendVBF("left", histos_1D[plot_index][fs][cat][Settings::AllData],
//                                        histos_1D[plot_index][fs][cat][Settings::H125VBF],
//                                        histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VH + ttH
//                                        histos_1D[plot_index][fs][cat][Settings::qqZZ],
//                                        histos_1D[plot_index][fs][cat][Settings::ggZZ],
//                                        histos_1D_ZX[plot_index][fs][cat]);
//   }
//   else if ( plot_index == Settings::D2jet_M4L118130 || plot_index == Settings::D2jet )
//   {
//      legend  = CreateLegendVBF("right", histos_1D[plot_index][fs][cat][Settings::AllData],
//                                         histos_1D[plot_index][fs][cat][Settings::H125VBF],
//                                         histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VH + ttH
//                                         histos_1D[plot_index][fs][cat][Settings::qqZZ],
//                                         histos_1D[plot_index][fs][cat][Settings::ggZZ],
//                                         histos_1D_ZX[plot_index][fs][cat]);
//   }
//   else if ( is_DVH_ )
//   {
//      legend  = CreateLegendVH("right", histos_1D[plot_index][fs][cat][Settings::AllData],
//                                        histos_1D[plot_index][fs][cat][Settings::H125VH],
//                                        histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VBF + ttH
//                                        histos_1D[plot_index][fs][cat][Settings::qqZZ],
//                                        histos_1D[plot_index][fs][cat][Settings::ggZZ],
//                                        histos_1D_ZX[plot_index][fs][cat]);
//   }
//   else
//   {
//      legend = CreateLegend((plot_index == Settings::MZ1 || plot_index == Settings::MZ1_M4L118130 || plot_index == Settings::MZ2 || plot_index == Settings::KD_M4L118130) ? "left" : "right",
//                             histos_1D[plot_index][fs][cat][Settings::AllData],
//                             histos_1D[plot_index][fs][cat][Settings::H125],
//                             histos_1D[plot_index][fs][cat][Settings::qqZZ],
//                             histos_1D[plot_index][fs][cat][Settings::ggZZ],
//                             histos_1D_ZX[plot_index][fs][cat]);
//   }
//
//   legend->Draw();
//
////===========
//// PLOT TEXT
////===========
//
//   TPaveText *text;
//
//   if ( plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::KD_M4L118130 || plot_index == Settings::MZ1_M4L118130 )
//   {
//      text = CreateCutText("right top", "118 < m_{4#font[12]{l}} < 130 GeV");
//      text->Draw();
//   }
//   else if ( plot_index == Settings::D2jet_M4L118130 || plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DZH_M4L118130 ||
//             plot_index == Settings::DVH_M4L118130   || plot_index == Settings::MZ2_M4L118130)
//   {
//      text = CreateCutText("left top", "118 < m_{4#font[12]{l}} < 130 GeV");
//      text->Draw();
//   }
//
////=================
//// CMS TEXT & LUMI
////=================
//
//   CMS_lumi *lumi = new CMS_lumi;
//   lumi->set_lumi(c, _lumi);
//
//   // Draw X-axis log scale
//   if ( plot_index == Settings::M4lMain )
//   {
//      stack->GetXaxis()->SetNdivisions(10);
//      stack->GetXaxis()->SetLabelSize(0);
//      DrawLogX(c, cat, fs);
//   }
//
//   _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_final_state.at(fs) + "_" + _s_category.at(cat);
//   SavePlots(c, _out_file_name);
//}
//=========================================================================================================


//=======================================
//void Histograms::SavePlots( TCanvas *c, TString name)
//{
//   c->SaveAs(name + ".pdf");
//   c->SaveAs(name + ".root");
////   c->SaveAs(name + ".C");
//   c->SaveAs(name + ".eps");
//   gSystem->Exec("convert -density 150 -quality 100 " + name + ".eps " + name + ".png");
//}
//=======================================







//============================================================================================================
//TLegend* Histograms::CreateLegend( string position, TH1F *data, TH1F *h125, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX )
//{
//   TLegend *leg;
//
//   if ( position == "right" )
//   {
//      leg = new TLegend(.71, .71, .91, .91);
//   }
//   else
//   {
//      leg = new TLegend(.21, .71, .41, .91);
//   }
//
//   leg->SetFillColor(0);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//
//   leg->AddEntry( data, "Data", "p E" );
//   leg->AddEntry( h125, "H(125)","f");
//   leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
//   leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
//   leg->AddEntry( ZX,   "Z+X", "f" );
//
//   return leg;
//}
//============================================================================================================


