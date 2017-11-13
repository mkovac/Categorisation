#ifndef Histograms_h
#define Histograms_h

// C++
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

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


   void SaveHistograms( TString );
   
//   void FillM4lZX( float, float, int, int );
//
//   void FillMZ1( float, float, float, int, int, int );
//   void FillMZ1ZX( float, float, float, int, int );
//
//   void FillMZ2( float, float, float, int, int, int );
//   void FillMZ2ZX( float, float, float, int, int );
//
//   void FillKD( float, float, float, int, int, int );
//   void FillKDZX( float, float, float, int, int );
//
//   void FillD1jet( float, float, float, int, int, int );
//   void FillD1jetZX( float, float, float, int, int );
//
//   void FillD2jet( float, float, float, int, int, int );
//   void FillD2jetZX( float, float, float, int, int );
//
//   void FillDWH( float, float, float, int, int, int );
//   void FillDWHZX( float, float, float, int, int );
//
//   void FillDZH( float, float, float, int, int, int );
//   void FillDZHZX( float, float, float, int, int );
//
//   void FillDVH( float, float, float, int, int, int );
//   void FillDVHZX( float, float, float, int, int );
//
//   void FillMZ1vsMZ2( float, float, float, float, int, int, int );
//
//   void FillVectors( float, float, float, int, float, float, float, float, float, int, int) ;
//   void FillDvsM4l( float, float, int, float, float, float, float, float, float, int, int, int );
//
//   void FillYields( float, float, int, int, int );
//
//   void SaveHistos( string );
//   void SaveYieldHistos( string );
//
//   void DeleteHistos();
//   void DeleteYieldsHistos();
//
//   void FillInclusive();
//   void FillInclusiveYields();
//
//   void SmoothHistograms();
//   void RenormalizeZX( vector<vector<float>> );
//
//   void GetHistos( TString );
//   void GetYieldsHistos( TString );
//
//   void plot_1D_single( TString, TString, TString, int, int );
//   void plot_1D_all_cat( TString, TString, TString );
//   void plot_1D_all_fs( TString, TString, TString );
//   void plot_2D_single( TString, TString, TString, int );
//   void plot_2D_error_single( TString, TString, TString, int );
//   void plot_2D_error_all_cat( TString , TString , TString );
//
//   void FillYieldGraphs( float, float );
//   void PrepareYamlFiles( TString , float , float, vector<vector<float>> );
//   void PrintYields( vector<vector<float>> );
//   void PrintYields( float, float, vector<vector<float>> );
//   void PrintLatexTables( float, float, vector<vector<float>> );
//   void setColZGradient_OneColor( int , bool );
//   void MakeZXShape( vector<vector<float>>, int );
//   void DrawLogX( TCanvas *, int, int );
//   void MakeCOLZGrey( bool );
//   void SavePlots( TCanvas *, TString );
//   void Rebin( THStack * );
//   void ChangeYaxisTitle( THStack * );
//   int SetProcess( int, int );
//   int SetPlotName( TString );
//   float SetMassPoint( int );
//
//   bool GetVarLogX( TString );
//   bool GetVarLogY( TString );
//
//   TLegend *CreateLegend( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* );
//   TLegend *CreateLegendVBF( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F* );
//   TLegend *CreateLegendVH( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F* );
//   TLegend *CreateLegendttH( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F* );
//   TLegend *Create2DLegend( string, TH2F*, TH2F*, TH2F* );
//   TLegend *Create2DErrorLegend( string, TGraphErrors*, TGraphErrors*, TGraphErrors* );
//   TLegend *Create2DLegendAllCat( string, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors* );
//
//   TPaveText *CreateCutText( string, TString);
//   TPaveText *CreateCatText( string, TString);
//
//   TLine *CreateDashedLine( int );
   

private:
      
   float lumi_;
   TString histo_name_, histo_label_;
   vector<TString> s_process_, s_sort_, s_gen_ch_;
   
//==============
// 1D histograms
//==============
   TH1F *h_pt_gen_H_lep_in_eta_acc_[Counters::num_of_processes][Counters::num_of_gen_channels][4];
   TH1F *h_pt_reco_bc_in_sig_reg_and_pass_triger_[Counters::num_of_processes][Counters::num_of_gen_channels][4];
   TH1F *h_eta_gen_H_lep_in_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels][4];
   TH1F *h_eta_reco_bc_in_sig_reg_and_pass_triger_[Counters::num_of_processes][Counters::num_of_gen_channels][4];

   TH1F *h_num_reco_H_ele_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   TH1F *h_num_reco_H_mu_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   TH1F *h_num_reco_H_lep_in_eta_pt_acc_[Counters::num_of_processes][Counters::num_of_gen_channels];
   
   TH1F *h_gen_H_pt_[Counters::num_of_processes];
   TH1F *h_gen_H_eta_[Counters::num_of_processes];
   TH2F *h_gen_H_eta_vs_pt_[Counters::num_of_processes];

   
};
#endif
