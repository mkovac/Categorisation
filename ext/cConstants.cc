#include "../include/cConstants.h"

#include <cassert>
#include <cmath>
#include <iostream>

cConstantSpline::cConstantSpline(const TString& filename) : filename_(filename), f_(nullptr), spline_(0) {}

void cConstants(){}

void cConstantSpline::initspline() {
   if (!spline_) {
      f_.reset(TFile::Open(filename_));
      spline_ = (TSpline3*)f_->Get("sp_gr_varTrue_Constant_Smooth");
   }
   assert(spline_);
}

double cConstantSpline::eval(double ZZMass) {
   initspline();
   return spline_->Eval(ZZMass);
}

cConstantSpline DbkgkinSpline2e2mu("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_Dbkgkin_2e2mu13TeV.root");
cConstantSpline DbkgkinSpline4e("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_Dbkgkin_4e13TeV.root");
cConstantSpline DbkgkinSpline4mu("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_Dbkgkin_4mu13TeV.root");
cConstantSpline DVBF2jetsSpline("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
cConstantSpline DVBF1jetSpline("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_DjVBF13TeV.root");
cConstantSpline DZHhSpline("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_DjjZH13TeV.root");
cConstantSpline DWHhSpline("/afs/cern.ch/work/m/mkovac/CMS/RUN_2/Data/Moriond_2017/c_constants/SmoothKDConstant_m4l_DjjWH13TeV.root");

extern "C" float getDVBF2jetsConstant(float ZZMass){
   return DVBF2jetsSpline.eval(ZZMass);
}
extern "C" float getDVBF1jetConstant(float ZZMass){
   return DVBF1jetSpline.eval(ZZMass);
}
extern "C" float getDWHhConstant(float ZZMass){
   return DWHhSpline.eval(ZZMass);
}
extern "C" float getDZHhConstant(float ZZMass){
   return DZHhSpline.eval(ZZMass);
}

extern "C" float getDVBF2jetsWP(float ZZMass, bool useQGTagging){
   if (useQGTagging) {
      assert(0);
      return 0.363;
   }
   else
      return 0.475598;
}
extern "C" float getDVBF1jetWP(float ZZMass, bool useQGTagging){
   if (useQGTagging) {
      assert(0);
      return 0.716;
   }
   else
      return 0.374228;
}
extern "C" float getDWHhWP(float ZZMass, bool useQGTagging){
   if (useQGTagging) {
      assert(0);
      return 0.965;
   }
   else
      //return 0.961587;
      return 0.941161;
}
extern "C" float getDZHhWP(float ZZMass, bool useQGTagging){
   if (useQGTagging) {
      assert(0);
      return 0.9952;
   }
   else
      //return 0.942319;
      return 0.941161;
}

extern "C" float getDVBF2jetsConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP) {
   float oldc = getDVBF2jetsConstant(ZZMass);
   float oldWP = getDVBF2jetsWP(ZZMass, useQGTagging);
   return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP));
}
extern "C" float getDVBF1jetConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP) {
   float oldc = getDVBF1jetConstant(ZZMass);
   float oldWP = getDVBF1jetWP(ZZMass, useQGTagging);
   return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP));
}

extern "C" float getDWHhConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP) {
   float oldc = getDWHhConstant(ZZMass);
   float oldWP = getDWHhWP(ZZMass, useQGTagging);
   return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP));
}

extern "C" float getDZHhConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP) {
   float oldc = getDZHhConstant(ZZMass);
   float oldWP = getDZHhWP(ZZMass, useQGTagging);
   return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP));
}


extern "C" float getDbkgkinConstant(int ZZflav, float ZZMass){ // ZZflav==id1*id2*id3*id4
   if (abs(ZZflav)==11*11*11*11 || abs(ZZflav)==2*11*11*11*11 || abs(ZZflav)==2*11*11*2*11*11) return DbkgkinSpline4e.eval(ZZMass);
   if (abs(ZZflav)==11*11*13*13 || abs(ZZflav)==2*11*11*13*13 || abs(ZZflav)==2*11*11*2*13*13) return DbkgkinSpline2e2mu.eval(ZZMass);
   if (abs(ZZflav)==13*13*13*13 || abs(ZZflav)==2*13*13*13*13 || abs(ZZflav)==2*13*13*2*13*13) return DbkgkinSpline4mu.eval(ZZMass);
   std::cout << "Invalid ZZflav " << ZZflav << std::endl; assert(0); return 0;
}
extern "C" float getDbkgConstant(int ZZflav, float ZZMass){
   return getDbkgkinConstant(ZZflav, ZZMass);
}
