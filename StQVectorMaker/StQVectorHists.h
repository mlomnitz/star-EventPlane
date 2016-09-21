#ifndef __qVectorHists__
#define __qVectorHists__

#include "TString.h"

class TH2F;
class TH3F;
class TProfile;
class TFile; 

class StQVectorHists {
 public:
  StQVectorHists(int harmonic);
  ~StQVectorHists();
  void addTrack(int mCent, float eta, float qx, float qy);
  void addEventPlane(int mCent, float Qx, float Qy);
  void Write(TFile *oFile);
 private:
  int mHarmonic;
  //
  TProfile*  prfQxCentEtaPlus;
  TProfile*  prfQyCentEtaPlus;
  TProfile*  prfQxCentEtaMinus;
  TProfile*  prfQyCentEtaMinus;
  
  TH2F*      hEventPlaneCent;
  TH3F*      hQyQxCent;

  ClassDef(StQVectorHists,0);
};
#endif
