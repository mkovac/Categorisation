// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include "Categorisation.h"

using namespace std;

int main( int argc, char *argv[] )
{
   gROOT->ProcessLine(".L ./ext/setTDRStyle_cpp.so");
   gROOT->ProcessLine("setTDRStyle();"); // Needed already here to save histos with markers
   
   TString path            = "/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/";
   TString file_name       = "/ZZ4lAnalysis.root";
   
   TString Data        = path + "AllData_2017" + file_name;
   TString ggH125      = path + "ggH125"       + file_name;
   TString VBFH125     = path + "VBFH125"      + file_name;
   TString WpH125      = path + "WplusH125"    + file_name;
   TString WmH125      = path + "WminusH125"   + file_name;
   TString ZH125       = path + "ZH125"        + file_name;
   TString bbH125      = path + "bbH125"       + file_name;
   TString ttH125      = path + "ttH125"       + file_name;
   TString ZZTo4l      = path + "ZZTo4l"       + file_name;
   TString ggZZ4e      = path + "ggTo4e"       + file_name;
   TString ggZZ4mu     = path + "ggTo4mu"      + file_name;
   TString ggZZ4tau    = path + "ggTo4tau"     + file_name;
   TString ggZZ2e2mu   = path + "ggTo2e2mu"    + file_name;
   TString ggZZ2e2tau  = path + "ggTo2e2tau"   + file_name;
   TString ggZZ2mu2tau = path + "ggTo2mu2tau"  + file_name;

   Categorisation *cat = new Categorisation(35.86706);

   cat->MakeHistograms(ggH125);
   
   cat->SaveHistograms("Histograms.root");

   
   delete cat;
}
