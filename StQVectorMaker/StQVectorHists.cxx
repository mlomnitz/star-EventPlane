#include <iostream>

#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TVector2.h" 
#include "TFile.h"
#include "TMath.h"

#include "StQVectorHists.h"

ClassImp(StQVectorHists);

StQVectorHists::StQVectorHists(int harmonic) : mHarmonic(harmonic)
{
  float PI = TMath::Pi();
  
  prfQxCentEtaPlus = new TProfile(Form("prfQxCentEtaPlus_v%i",harmonic),Form("prfQxCentEtaPlus_v%i",harmonic),9,0,9);
  prfQyCentEtaPlus = new TProfile(Form("prfQyCentEtaPlus_v%i",harmonic),Form("prfQyCentEtaPlus_v%i",harmonic),9,0,9);
  prfQxCentEtaMinus = new TProfile(Form("prfQxCentEtaMinus_v%i",harmonic),Form("prfQxCentEtaMinus_v%i",harmonic),9,0,9);
  prfQyCentEtaMinus = new TProfile(Form("prfQyCentEtaMinus_v%i",harmonic),Form("prfQyCentEtaMinus_v%i",harmonic),9,0,9);
  
  hEventPlaneCent = new TH2F(Form("hEventPlaneCent_v%i",harmonic),Form("hEventPlaneCent_v%i",harmonic),9,0,9,60,0,2.0*PI/harmonic);
  hQyQxCent = new TH3F(Form("hQyQxCent_v%i",harmonic),Form("hQyQxCent_v%i",harmonic), 9,0,9,1000,-50,50,1000,-50,50);
}
StQVectorHists::~StQVectorHists()
{ } 
void StQVectorHists::addTrack(int mCent, float eta, float qx, float qy)
{
  if(eta>0){
    prfQxCentEtaPlus->Fill(mCent, qx);
    prfQyCentEtaPlus->Fill(mCent, qy);
  }
  else{
    prfQxCentEtaMinus->Fill(mCent, qx);
    prfQyCentEtaMinus->Fill(mCent, qy);
  }
}
void StQVectorHists::addEventPlane(int mCent, float Qx, float Qy)
{
  TVector2 Q(Qx,Qy);
  float eventPlane = Q.Phi()/mHarmonic;
  hQyQxCent->Fill(mCent, Qx, Qy);
  if(Q.Mod()>0)
    hEventPlaneCent->Fill(mCent, eventPlane);
}
void StQVectorHists::Write(TFile* oFile)
{
  if( !oFile)
    std::cout<<"Warning, no output file"<<std::endl;
  oFile->cd();
  prfQxCentEtaPlus->Write();
  prfQyCentEtaPlus->Write();
  prfQxCentEtaMinus->Write();
  prfQyCentEtaMinus->Write();
  hQyQxCent->Write();
  hEventPlaneCent->Write();
}
