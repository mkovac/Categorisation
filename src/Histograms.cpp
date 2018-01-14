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
   
   // Associated decay
   s_assoc_decay_.push_back("WH_4_0");
   s_assoc_decay_.push_back("WH_4_1");
   s_assoc_decay_.push_back("ZH_4_0");
   s_assoc_decay_.push_back("ZH_4_2");
   s_assoc_decay_.push_back("ZH_2_2");
   s_assoc_decay_.push_back("bbH_0");
   s_assoc_decay_.push_back("bbH_1");
   s_assoc_decay_.push_back("ttH_4_0");
   s_assoc_decay_.push_back("ttH_4_1");
   s_assoc_decay_.push_back("ttH_4_2");
   s_assoc_decay_.push_back("ttH_2_2");

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
   
   
   // ROC names
   s_ROC_.push_back("DiJetFisher_qqH_ggH");
   s_ROC_.push_back("DiJetFisher_qqH_qqZZ");
   s_ROC_.push_back("Pvbf_qqH_ggH");
   s_ROC_.push_back("Pvbf_qqH_qqZZ");
   s_ROC_.push_back("Phjj_qqH_ggH");
   s_ROC_.push_back("Phjj_qqH_qqZZ");
   s_ROC_.push_back("Phjj_WH_ggH");
   s_ROC_.push_back("Phjj_ZH_ggH");
   s_ROC_.push_back("Pvbf1j_qqH_ggH");
   s_ROC_.push_back("Phj_qqH_ggH");
   s_ROC_.push_back("Pwhhadr_WH_ggH");
   s_ROC_.push_back("Pzhhadr_ZH_ggH");
   s_ROC_.push_back("D2jVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D1jVbfHj_qqH_ggH");
   s_ROC_.push_back("D2jWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("Pqj1Pqj2_qqH_ggH");
   s_ROC_.push_back("Pgj1Pgj2_qqH_ggH");
   s_ROC_.push_back("D2jqg_qqH_ggH");
   s_ROC_.push_back("Dqgj1Dqgj2_qqH_ggH");
   s_ROC_.push_back("Pqj1Pqj2_qqH_qqZZ");
   s_ROC_.push_back("Pgj1Pgj2_qqH_qqZZ");
   s_ROC_.push_back("D2jqg_qqH_qqZZ");
   s_ROC_.push_back("Dqgj1Dqgj2_qqH_qqZZ");
   s_ROC_.push_back("Pq_qqH_ggH");
   s_ROC_.push_back("Pg_qqH_ggH");
   s_ROC_.push_back("D1jqg_qqH_ggH");
   s_ROC_.push_back("Pqj1Pqj2_WH_ggH");
   s_ROC_.push_back("Pgj1Pgj2_WH_ggH");
   s_ROC_.push_back("D2jqg_WH_ggH");
   s_ROC_.push_back("Dqgj1Dqgj2_WH_ggH");
   s_ROC_.push_back("Pqj1Pqj2_ZH_ggH");
   s_ROC_.push_back("Pgj1Pgj2_ZH_ggH");
   s_ROC_.push_back("D2jqg_ZH_ggH");
   s_ROC_.push_back("Dqgj1Dqgj2_ZH_ggH");
   s_ROC_.push_back("D2jMelaQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaD2jQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaD2jQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D1jMelaQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D1jMelaD1jQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D2jMelaQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaD2jQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaQGZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("D2jMelaD2jQGZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("RatioPvbfPhjj_qqH_ggH");
   s_ROC_.push_back("RatioPqj1Pqj2Pgj1Pgj2_qqH_ggH");
   s_ROC_.push_back("RatioPvbfPqj1Pqj2PhjjPgj1Pgj2_qqH_ggH");
   s_ROC_.push_back("D2jMelaExpQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaSqQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaSqrtQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaCbrtQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaQrrtQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaQnrtQGVbfHjj_qqH_ggH");
   s_ROC_.push_back("D2jMelaExpQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaSqQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaSqrtQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaCbrtQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaQrrtQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D2jMelaQnrtQGVbfHjj_qqH_qqZZ");
   s_ROC_.push_back("D1jMelaSqrtQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D1jMelaCbrtQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D1jMelaQrrtQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D1jMelaQnrtQGVbfHj_qqH_ggH");
   s_ROC_.push_back("D2jMelaSqrtQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaCbrtQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaQrrtQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaQnrtQGWHHadrHjj_WH_ggH");
   s_ROC_.push_back("D2jMelaSqrtQGZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("D2jMelaCbrtQGZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("D2jMelaQrrtQGZHHadrHjj_ZH_ggH");
   s_ROC_.push_back("D2jMelaQnrtQGZHHadrHjj_ZH_ggH");
   
   
   //Prepare variable pairs
   var_pair.PrepareVarPair("M4l_vs_D_kin_bkg","m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", 100, 50, 850, 20, 0, 1);
   var_pair.PrepareVarPair("MZ2_vs_D_kin_bkg", "m_{Z_{2}} (GeV)", "D_{bkg}^{kin}", 75, 0, 150, 20, 0, 1);
   var_pair.PrepareVarPair("D2jVbfHjj_vs_D2jqg", "D_{VBF-2j}^{ME}", "D_{2jets}^{q/g}", 51, 0, 1.02, 25, -0.04, 0.96);


   //Prepare variables: n_bins, x_min, x_max, x_max_roc
   var.PrepareVar("M4l", "m_{4#font[12]{l}} (GeV)", 100, 50, 85, 850);
   var.PrepareVar("M4l2", "m_{4#font[12]{l}} (GeV)", 70, 105, 140, 140);
   var.PrepareVar("MZ1", "m_{Z_{1}} (GeV)", 75, 0, 150, 150);
   var.PrepareVar("MZ2", "m_{Z_{2}} (GeV)", 75, 0, 150, 150);
   var.PrepareVar("Dkinbkg", "D_{bkg}^{kin}", 20, 0, 1, 1);
   var.PrepareVar("DjetFisher", "D^{Fisher}", 50, 0, 2, 20);
   var.PrepareVar("Pvbf", "P_{VBF}^{MELA}", 50, -2, 98, 10000);
   var.PrepareVar("Phjj", "P_{ggH+2j}^{MELA}", 50, -5, 245, 10000);
   var.PrepareVar("Pvbf1j", "P_{VBF 1j}^{MELA}", 50, -10, 490, 10000);
   var.PrepareVar("Phj", "P_{ggH+1j}^{MELA}", 50, -10, 490, 10000);
   var.PrepareVar("Pwhhadr", "P_{WH-h}^{MELA}", 50, -10, 490, 1000000);
   var.PrepareVar("Pzhhadr", "P_{ZH-h}^{MELA}", 50, -10, 490, 1000000);
   var.PrepareVar("Pwhlept", "P_{WH-l}^{MELA}", 50, -1, 49, 1000);
   var.PrepareVar("Pzhlept", "P_{ZH-l}^{MELA}", 50, -1, 49, 1000);
   var.PrepareVar("D2jVbfHjj", "D_{VBF-2j}^{ME}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jVbfHj", "D_{VBF-1j}^{ME}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jWHHadrHjj", "D_{WH-hadr.}^{ME}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jZHHadrHjj", "D_{ZH-hadr.}^{ME}", 51, 0, 1.02, 1.02);
   var.PrepareVar("Pqj1", "P_{q}(j_{1})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pgj1", "P_{g}(j_{1})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pqj1Pqj2", "P_{q}(j_{1})*P_{q}(j_{2})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pgj1Pgj2", "P_{g}(j_{1})*P_{g}(j_{2})", 25, -0.04, 0.96, 100);
   var.PrepareVar("D2jqg", "D_{2jets}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("Dqgj1Dqgj2", "D_{1jet}^{q/g}(j_{1})*D_{1jet}^{q/g}(j_{2})", 51, 0, 1.02, 1.02);
   var.PrepareVar("Pqj1VbfTopo", "P_{q}(j_{1})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pgj1VbfTopo", "P_{g}(j_{1})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pqj1Pqj2VbfTopo", "P_{q}(j_{1})*P_{q}(j_{2})", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pgj1Pgj2VbfTopo", "P_{g}(j_{1})*P_{g}(j_{2})", 25, -0.04, 0.96, 100);
   var.PrepareVar("D2jqgVbfTopo", "D_{2jets}^{q/g}", 26, 0, 1.04, 1.04);
   var.PrepareVar("Pq", "P_{q}", 25, -0.04, 0.96, 100);
   var.PrepareVar("Pg", "P_{g}", 25, -0.04, 0.96, 100);
   var.PrepareVar("D1jqg", "D_{1jet}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQGVbfHjj", "D_{VBF-2j}^{ME q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaD2jQGVbfHjj", "D_{VBF-2j}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaQGVbfHj", "D_{VBF-1j}^{ME q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaD1jQGVbfHj", "D_{VBF-1j}^{ME}*D_{1jet}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQGWHHadrHjj", "D_{WH-hadr.}^{ME q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaD2jQGWHHadrHjj", "D_{WH-hadr.}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQGZHHadrHjj", "D_{ZH-hadr.}^{ME q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaD2jQGZHHadrHjj", "D_{ZH-hadr.}^{ME}*D_{2jets}^{q/g}", 51, 0, 1.02, 1.02);
   var.PrepareVar("RatioPvbfPhjj", "P_{VBF}^{MELA} / P_{ggH+2j}^{MELA}", 50, 0, 5, 100);
   var.PrepareVar("RatioPqj1Pqj2Pgj1Pgj2", "[P_{q}(j_{1})*P_{q}(j_{2})] / [P_{g}(j_{1})*P_{g}(j_{2})]", 50, 0, 5, 100);
   var.PrepareVar("RatioPvbfPqj1Pqj2PhjjPgj1Pgj2", "[P_{VBF}^{MELA}*P_{q}(j_{1})*P_{q}(j_{2})]/[P_{ggH+2j}^{MELA}*P_{g}(j_{1})*P_{g}(j_{2})]", 50, 0, 5, 100);
   var.PrepareVar("D2jMelaExpQGVbfHjj", "D_{VBF-2j}^{ME exp(q/g)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaSqQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^2}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaSqrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/2)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaCbrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/3)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQrrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/4)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQnrtQGVbfHjj", "D_{VBF-2j}^{ME (q/g)^(1/5)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaSqrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/2)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaCbrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/3)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaQrrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/4)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D1jMelaQnrtQGVbfHj", "D_{VBF-1j}^{ME (q/g)^(1/5)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaSqrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/2)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaCbrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/3)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQrrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/4)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQnrtQGWHHadrHjj", "D_{WH-hadr.}^{ME (q/g)^(1/5)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaSqrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/2)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaCbrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/3)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQrrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/4)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("D2jMelaQnrtQGZHHadrHjj", "D_{ZH-hadr.}^{ME (q/g)^(1/5)}", 51, 0, 1.02, 1.02);
   var.PrepareVar("Pt4l", "p_{T}^{4l} (GeV)", 50, 0, 500, 500);
   var.PrepareVar("NGenLep", "# gen leptons", 7, 0, 7, 7);
   var.PrepareVar("NGenLepInEtaPtAcc", "# gen leptons in acceptance", 7, 0, 7, 7);
   var.PrepareVar("NGenLepNotInEtaPtAcc", "# gen leptons not in acceptance", 7, 0, 7, 7);
   var.PrepareVar("NGenHLepNotInEtaPtAcc", "# gen leptons from H not in acceptance", 5, 0, 5, 5);
   var.PrepareVar("NGenAssocLepNotInEtaPtAcc", "# gen associated leptons not in acceptance", 3, 0, 3, 3);
   var.PrepareVar("NGenLepMinusNGoodLep", "# gen leptons - # good leptons", 6, -3, 3, 3);
   var.PrepareVar("NGenLepInEtaPtAccMinusNGoodLep", "# gen leptons in acceptance - # good leptons", 6, -3, 3, 3);
   var.PrepareVar("NExtraLep", "number of additional leptons", 6, 0, 6, 6);
   var.PrepareVar("NExtraZ", "number of additional #font[12]{l}^{+}#font[12]{l}^{-} pairs", 4, 0, 4, 4);
   var.PrepareVar("NJets", "number of selected jets", 15, 0, 15, 15);
   var.PrepareVar("NBtaggedJets", "number of selected b-tagged jets", 6, 0, 6, 6);
   var.PrepareVar("MET", "E_{T}^{miss.}", 100, 0, 200, 200);
   
   //Link ROCs and variables
   ROC_to_variables.push_back(Counters::DiJetFisher);
   ROC_to_variables.push_back(Counters::DiJetFisher);
   ROC_to_variables.push_back(Counters::Pvbf);
   ROC_to_variables.push_back(Counters::Pvbf);
   ROC_to_variables.push_back(Counters::Phjj);
   ROC_to_variables.push_back(Counters::Phjj);
   ROC_to_variables.push_back(Counters::Phjj);
   ROC_to_variables.push_back(Counters::Phjj);
   ROC_to_variables.push_back(Counters::Pvbf1j);
   ROC_to_variables.push_back(Counters::Phj);
   ROC_to_variables.push_back(Counters::Pwhhadr);
   ROC_to_variables.push_back(Counters::Pzhhadr);
   ROC_to_variables.push_back(Counters::D2jVbfHjj);
   ROC_to_variables.push_back(Counters::D2jVbfHjj);
   ROC_to_variables.push_back(Counters::D1jVbfHj);
   ROC_to_variables.push_back(Counters::D2jWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jZHHadrHjj);
   ROC_to_variables.push_back(Counters::Pqj1Pqj2);
   ROC_to_variables.push_back(Counters::Pgj1Pgj2);
   ROC_to_variables.push_back(Counters::D2jqg);
   ROC_to_variables.push_back(Counters::Dqgj1Dqgj2);
   ROC_to_variables.push_back(Counters::Pqj1Pqj2);
   ROC_to_variables.push_back(Counters::Pgj1Pgj2);
   ROC_to_variables.push_back(Counters::D2jqg);
   ROC_to_variables.push_back(Counters::Dqgj1Dqgj2);
   ROC_to_variables.push_back(Counters::Pq);
   ROC_to_variables.push_back(Counters::Pg);
   ROC_to_variables.push_back(Counters::D1jqg);
   ROC_to_variables.push_back(Counters::Pqj1Pqj2);
   ROC_to_variables.push_back(Counters::Pgj1Pgj2);
   ROC_to_variables.push_back(Counters::D2jqg);
   ROC_to_variables.push_back(Counters::Dqgj1Dqgj2);
   ROC_to_variables.push_back(Counters::Pqj1Pqj2);
   ROC_to_variables.push_back(Counters::Pgj1Pgj2);
   ROC_to_variables.push_back(Counters::D2jqg);
   ROC_to_variables.push_back(Counters::Dqgj1Dqgj2);
   ROC_to_variables.push_back(Counters::D2jMelaQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaD2jQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaD2jQGVbfHjj);
   ROC_to_variables.push_back(Counters::D1jMelaQGVbfHj);
   ROC_to_variables.push_back(Counters::D1jMelaD1jQGVbfHj);
   ROC_to_variables.push_back(Counters::D2jMelaQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaD2jQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQGZHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaD2jQGZHHadrHjj);
   ROC_to_variables.push_back(Counters::RatioPvbfPhjj);
   ROC_to_variables.push_back(Counters::RatioPqj1Pqj2Pgj1Pgj2);
   ROC_to_variables.push_back(Counters::RatioPvbfPqj1Pqj2PhjjPgj1Pgj2);
   ROC_to_variables.push_back(Counters::D2jMelaExpQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaSqQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaSqrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaCbrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQrrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQnrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaExpQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaSqQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaSqrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaCbrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQrrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQnrtQGVbfHjj);
   ROC_to_variables.push_back(Counters::D1jMelaSqrtQGVbfHj);
   ROC_to_variables.push_back(Counters::D1jMelaCbrtQGVbfHj);
   ROC_to_variables.push_back(Counters::D1jMelaQrrtQGVbfHj);
   ROC_to_variables.push_back(Counters::D1jMelaQnrtQGVbfHj);
   ROC_to_variables.push_back(Counters::D2jMelaSqrtQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaCbrtQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQrrtQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQnrtQGWHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaSqrtQGZHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaCbrtQGZHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQrrtQGZHHadrHjj);
   ROC_to_variables.push_back(Counters::D2jMelaQnrtQGZHHadrHjj);
   

   
   
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
         histo_label_ = ";;# of expected events";
         histo_name_ = "h_BC_in_SR_basket_" +  s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
         h_BC_in_SR_basket_[i_proc][i_rc] = new TH1F(histo_name_, histo_label_, 8, 0, 8 );
         
         
         for ( int i_var = 0; i_var < Counters::num_of_vars; i_var++ )
         {
            histo_label_ = ";" + var.vec_var[i_var].label + ";" + "# of events";

            histo_name_ = "h_bc_in_sig_reg_" + var.vec_var[i_var].name + "_" +  s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
            h_bc_in_sig_reg_[i_var][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                        var.vec_var[i_var].x_min,
                                                                                        var.vec_var[i_var].x_max);
            
            for ( int i_ms = 0; i_ms < Counters::num_of_H_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_H_leps_" + var.vec_var[i_var].name + "_" + s_H_lep_match_status_.at(i_ms) + "_"
                                                             + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_H_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                                              var.vec_var[i_var].x_min,
                                                                                                              var.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_all_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_all_leps_" + var.vec_var[i_var].name + "_" + s_all_lep_match_status_.at(i_ms) + "_"
                                                               + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_all_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                                                var.vec_var[i_var].x_min,
                                                                                                                var.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_WH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_WH_leps_" + var.vec_var[i_var].name + "_" + s_WH_lep_match_status_.at(i_ms) + "_"
                                                              + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_WH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                                               var.vec_var[i_var].x_min,
                                                                                                               var.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_ZH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_ZH_leps_" + var.vec_var[i_var].name + "_" + s_ZH_lep_match_status_.at(i_ms) + "_"
                                                              + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_ZH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                                               var.vec_var[i_var].x_min,
                                                                                                               var.vec_var[i_var].x_max);
            }
            
            for ( int i_ms = 0; i_ms < Counters::num_of_ttH_lep_match_statuses; i_ms++ )
            {
               histo_name_ = "h_bc_in_sig_reg_match_ttH_leps_" + var.vec_var[i_var].name + "_" + s_ttH_lep_match_status_.at(i_ms) + "_"
                                                               + s_process_.at(i_proc) + "_" + s_reco_ch_.at(i_rc);
               h_bc_in_sig_reg_match_ttH_leps_[i_var][i_ms][i_proc][i_rc] = new TH1F(histo_name_, histo_label_, var.vec_var[i_var].n_bins,
                                                                                                                var.vec_var[i_var].x_min,
                                                                                                                var.vec_var[i_var].x_max);
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


   for ( int i_rc = 0; i_rc < Counters::num_of_reco_ch; i_rc++ )
   {
      for ( int i_assoc = 0; i_assoc < Counters::num_of_associated_decays; i_assoc++ )
      {
         histo_label_ = ";;# of expected events";
         histo_name_ = "h_BC_in_SR_basket_assoc_dec_" +  s_assoc_decay_.at(i_assoc) + "_" + s_reco_ch_.at(i_rc);
         h_BC_in_SR_basket_assoc_dec_[i_assoc][i_rc] = new TH1F(histo_name_, histo_label_, 8, 0, 8 );
      }
   }
   
   
   for ( int i_roc = 0; i_roc < Counters::num_of_ROCs; i_roc++ )
   {
      histo_label_ = ";" + var.vec_var[ROC_to_variables[i_roc]].label + ";" + "# of events";
     
      histo_name_ = "h_ROC_sig_" + s_ROC_.at(i_roc);
      h_ROC_sig_[i_roc] = new TH1F(histo_name_, histo_label_, var.vec_var[ROC_to_variables[i_roc]].n_bins,
                                                              var.vec_var[ROC_to_variables[i_roc]].x_min,
                                                              var.vec_var[ROC_to_variables[i_roc]].x_max_roc);
      
      histo_name_ = "h_ROC_bkg_" + s_ROC_.at(i_roc);
      h_ROC_bkg_[i_roc] = new TH1F(histo_name_, histo_label_, var.vec_var[ROC_to_variables[i_roc]].n_bins,
                                                              var.vec_var[ROC_to_variables[i_roc]].x_min,
                                                              var.vec_var[ROC_to_variables[i_roc]].x_max_roc);
   } // end i_roc
   
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



//=======================================================================
void Histograms::FillSigROC( float var_value, int ROC_num, float weight )
{
   h_ROC_sig_[ROC_num]->Fill(var_value, weight);
}
//=======================================================================



//=======================================================================
void Histograms::FillBkgROC( float var_value, int ROC_num, float weight )
{
   h_ROC_bkg_[ROC_num]->Fill(var_value, weight);
}
//=======================================================================



//================================================================================================
void Histograms::FillBasketHistograms( int var_value, float weight, int proc, int rc_1, int rc_2 )
{
   h_BC_in_SR_basket_[proc][rc_1]->Fill(var_value, weight);
   h_BC_in_SR_basket_[proc][rc_1]->Fill(7., weight);
   h_BC_in_SR_basket_[proc][rc_2]->Fill(var_value, weight);
   h_BC_in_SR_basket_[proc][rc_2]->Fill(7., weight);
}
//================================================================================================



//================================================================================================
void Histograms::FillBasketHistogramsAssoc( int var_value, float weight, int assoc_dec, int rc_1, int rc_2 )
{
   h_BC_in_SR_basket_assoc_dec_[assoc_dec][rc_1]->Fill(var_value, weight);
   h_BC_in_SR_basket_assoc_dec_[assoc_dec][rc_1]->Fill(7., weight);
   h_BC_in_SR_basket_assoc_dec_[assoc_dec][rc_2]->Fill(var_value, weight);
   h_BC_in_SR_basket_assoc_dec_[assoc_dec][rc_2]->Fill(7., weight);
}
//================================================================================================



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
         h_BC_in_SR_basket_[i_proc][i_rc]->Write();
         
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
   
   
   for ( int i_rc = 0; i_rc < Counters::num_of_reco_ch; i_rc++ )
   {
      for ( int i_assoc = 0; i_assoc < Counters::num_of_associated_decays; i_assoc++ )
      {
         h_BC_in_SR_basket_assoc_dec_[i_assoc][i_rc]->Write();
      }
   }
   
   
   for ( int i_roc = 0; i_roc < Counters::num_of_ROCs; i_roc++ )
   {
      h_ROC_sig_[i_roc]->Write();
      h_ROC_bkg_[i_roc]->Write();
   } // end i_roc
   
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



//=============================================
//void Histograms::PlotCategories()
//{
//   TCanvas *c_BC_in_SR_baskets_all = new TCanvas("c_BC_in_SR_baskets_all", "c_BC_in_SR_baskets_all", 1500, 500 );
//    DrawByProdmodes(cBCInSRBasketsAll,hBCInSRBaskets,-1);
//
//
//
//
//
//
//
//
//   
//
//
//    SaveCanvas(outDir,cBCInSRBasketsAll,tagOut);
//
//
//
//
//    TCanvas* cBCInSRBasketEfficiency[nProcesses];
//    TH1F* hProdModes[nProcesses]; for(int pr=0; pr<nProcesses; pr++) hProdModes[pr] = (TH1F*)hBCInSRBaskets[pr][FINALSTATE];
//    TH1F* hAssocDecays[nAssocDecays]; for(int a=0; a<nAssocDecays; a++) hAssocDecays[a] = (TH1F*)hBCInSRBasketsAssocDecays[a][FINALSTATE];
//    for(int pr=0; pr<nProcesses; pr++){
//      if(!isProcessed[pr]) continue;
//      cBCInSRBasketEfficiency[pr] = new TCanvas(Form("cBCInSR_basketEfficiency_%s",sProcess[pr].c_str()),Form("cBCInSR_basketEfficiency_%s",sProcess[pr].c_str()),500,longDim[BASKETLIST]);
//      DrawBasketEfficiencies(cBCInSRBasketEfficiency[pr],hProdModes[pr],pr,hAssocDecays,treatH2l2XAsBkgd?assocDecayName2:assocDecayName1,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumiText,false);
//      SaveCanvas(outDir,cBCInSRBasketEfficiency[pr],tagOut);
//    }
//
//    TCanvas* cBCInSRBasketPurity = new TCanvas("cBCInSR_basketPurity","cBCInSR_basketPurity",800,longDim[BASKETLIST]);
//    DrawBasketPurities(cBCInSRBasketPurity,hProdModes,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,false,treatH2l2XAsBkgd,basketLabel[BASKETLIST],labelMerge,lumiText);
//    SaveCanvas(outDir,cBCInSRBasketPurity,tagOut);
//    TCanvas* cBCInSRBasketPuritySplit = new TCanvas("cBCInSR_basketPuritySplit","cBCInSR_basketPuritySplit",800,longDim[BASKETLIST]);
//    DrawBasketPurities(cBCInSRBasketPuritySplit,hProdModes,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,true,treatH2l2XAsBkgd,basketLabel[BASKETLIST],labelMerge,lumiText);
//    SaveCanvas(outDir,cBCInSRBasketPuritySplit,tagOut);
//
//    if(isProcessed[qqZZ]){
//
//      TH1F* hSumSgnl = (TH1F*)hBCInSRBaskets[ggH][FINALSTATE]->Clone();
//      for(int pr=ggH; pr<nSignalProcesses; pr++){
//   if(!isProcessed[pr]) continue;
//   hSumSgnl->Add(hBCInSRBaskets[pr][FINALSTATE]);
//      }
//      TH1F* hSumBkgd = (TH1F*)hBCInSRBaskets[qqZZ][FINALSTATE]->Clone();
//      for(int pr=qqZZ; pr<nProcesses; pr++){
//   if(!isProcessed[pr]) continue;
//   hSumBkgd->Add(hBCInSRBaskets[pr][FINALSTATE]);
//      }
//      if(treatH2l2XAsBkgd){
//   hSumSgnl->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE],-1);
//   hSumSgnl->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE],-1);
//   hSumBkgd->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE]);
//   hSumBkgd->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE]);
//      }
//
//      TH1F* hDenom = (TH1F*)hSumSgnl->Clone();
//      hDenom->Add(hSumBkgd);
//      TH1F* hSOSPB = (TH1F*)hSumSgnl->Clone();
//      hSOSPB->Divide(hDenom);
//      TCanvas* cBCInSRBasketSOSPB = new TCanvas("cBCInSR_basketSOSPB","cBCInSR_basketSOSPB",500,longDim[BASKETLIST]);
//      DrawBasketSOSPB(cBCInSRBasketSOSPB,hSOSPB,basketLabel[BASKETLIST]);
//      SaveCanvas(outDir,cBCInSRBasketSOSPB,tagOut);
//
//    }
//
//    for(int pr=0; pr<nProcesses; pr++)
//      if(isProcessed[pr]) cout<<"process "<<sProcess[pr]<<": "<<0.5*hBCInSRBaskets[pr][FINALSTATE]->Integral()<<endl;
//
//    //*
//    for(int pr=0; pr<nProcesses; pr++)
//      for(int ba=0; ba<nBaskets; ba++)
//   if(isProcessed[pr]) cout<<" basket #"<<ba<<", process "<<sProcess[pr]<<": "<<hBCInSRBaskets[pr][FINALSTATE]->GetBinContent(ba+1)<<endl;
//    //*/
//
//  }
//
//
//}
//=============================================



//=========================================================================================================
//void DrawByProdMode(TCanvas *c, TH1F *histo[Counters::num_of_processes][Counters::num_of_reco_ch], int v)
//{
//
//  gStyle->SetOptTitle(0);
//
//  c->cd();
//  c->SetTicks(0, 0);
//
//  Float_t max = 0.;
//
//   for ( int i_proc = 0; i_proc < Counters::num_of_processes; i_proc++ )
//   {
//
//    histo[pr]->Scale(1/h[pr]->Integral(0,h[pr]->GetNbinsX()+1));
//    Float_t maxtemp = h[pr]->GetMaximum();
//    if(maxtemp>max) max = maxtemp;
//
//  }
//  Float_t cmax = logY ? 10.*max : (v==Dkinbkg?2.:v==D1jVbfHj?1.3:1.1)*max ;
//  Float_t cminlog = max / varMinFactor[v];
//
//  for(int pr=0; pr<nProcesses; pr++){
//    if(!isProcessed[pr]) continue;
//    if(SKIPPROCESSES && skipProcess[pr]) continue;
//
//    h[pr]->SetLineColor(processColor[pr]);
//    h[pr]->SetLineWidth(2);
//    h[pr]->SetFillStyle(0);
//    h[pr]->SetStats(0);
//    //if(!logY) h[pr]->SetMinimum(0);
//
//    if(pr==0){
//      h[pr]->SetMaximum(cmax);
//      h[pr]->SetMinimum(logY?cminlog:0.);
//      h[pr]->GetYaxis()->SetTitle("normalized to 1");
//      h[pr]->Draw("hist");
//    }else{
//      h[pr]->Draw("histsame");
//    }
//
//    lgd->AddEntry( h[pr], processLabel[pr].c_str(), "l" );
//  }
//
//  lgd->Draw();
//
//  if(v>=0 && varCutLabel[v]!="")
//    printInfo(varCutLabel[v],0.72,lgdLow-0.09,1.,lgdLow-0.04);
//
//  gPad->RedrawAxis();
//
//  writeExtraText = true;
//  extraText  = "Simulation";
//  lumi_sqrtS = "13 TeV";
//  CMS_lumi( c, 0, 0);
//
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


