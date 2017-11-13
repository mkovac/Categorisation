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
   
   s_sort_.push_back("first");
   s_sort_.push_back("second");
   s_sort_.push_back("third");
   s_sort_.push_back("fourth");
   
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

   
//   _s_final_state.push_back("4e");
//   _s_final_state.push_back("4mu");
//   _s_final_state.push_back("2e2mu");
//   _s_final_state.push_back("2mu2e");
//   _s_final_state.push_back("4l");
//
//   _s_category.push_back("UnTagged");
//   _s_category.push_back("VBF1jTagged");
//   _s_category.push_back("VBF2jTagged");
//   _s_category.push_back("VHLeptTagged");
//   _s_category.push_back("VHHadrTagged");
//   _s_category.push_back("ttHTagged");
//   _s_category.push_back("VHMETTagged");
//   _s_category.push_back("Inclusive");
//
//   _s_production_mode.push_back("ggH");
//   _s_production_mode.push_back("qqH");
//   _s_production_mode.push_back("WH");
//   _s_production_mode.push_back("ZH");
//   _s_production_mode.push_back("ttH");
   

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
   
   
      for ( int i_gc = 0; i_gc < Counters::num_of_gen_channels; i_gc++ )
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

         for ( int i_sort = 0; i_sort < 4; i_sort++ )
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
         } // end for
      } // end for
   } // end for
   
   
   
   
   
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
   
      for ( int i_gc = 0; i_gc < Counters::num_of_gen_channels; i_gc++ )
      {
         h_num_reco_H_ele_in_eta_pt_acc_[i_proc][i_gc]->Write();
         h_num_reco_H_mu_in_eta_pt_acc_[i_proc][i_gc]->Write();
         h_num_reco_H_lep_in_eta_pt_acc_[i_proc][i_gc]->Write();
      
         for ( int i_sort = 0; i_sort < 4; i_sort++ )
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


