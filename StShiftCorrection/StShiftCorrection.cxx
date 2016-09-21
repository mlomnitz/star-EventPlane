#define only__primary
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TProfile.h"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "PhysicalConstants.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "StShiftCorrection.h"

ClassImp(StShiftCorrection)

//-----------------------------------------------------------------------------
StShiftCorrection::StShiftCorrection(const char* name, StPicoDstMaker *picoMaker, StRefMultCorr* grefmultCorrUtil)
   : StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(NULL),  mPicoEvent(NULL), mgrefmultCorrUtil(grefmultCorrUtil),
  mAcceptEvent(false), mAcceptQvectorFile(false), mAcceptQvectorFiletmp(true),mCent(-1), mRunNumber(0), mBField(-999.),mVertexPos(-999, -999, -999)
{
  TFile* outFile = new TFile(Form("%s.qVectr_shift_corrected.root",name),"RECREATE");
  setFileOut(outFile);
  for( int ii = 0 ; ii<ShiftConstants::mNoEP; ++ii){
    prfQxCentEtaPlus[ii] = NULL;
    prfQyCentEtaPlus[ii] = NULL;
    prfQxCentEtaMinus[ii] = NULL;
    prfQxCentEtaMinus[ii] = NULL;
    for(int jj = 0 ; jj < ShiftConstants::mShiftMaxHarmonic; ++jj){
      prfShiftQxCent[ii][jj] = new TProfile(Form("prfShiftQxCent_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
					    Form("prfShiftQxCent_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
      prfShiftQyCent[ii][jj] = new TProfile(Form("prfShiftQyCent_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
					    Form("prfShiftQyCent_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
      //
      prfShiftQxCentEtaPlus[ii][jj] = new TProfile(Form("prfShiftQxCentEtaPlus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
						   Form("prfShiftQxCentEtaPlus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
      prfShiftQyCentEtaPlus[ii][jj] = new TProfile(Form("prfShiftQyCentEtaPlus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
						   Form("prfShiftQyCentEtaPlus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
      prfShiftQxCentEtaMinus[ii][jj] = new TProfile(Form("prfShiftQxCentEtaMinus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
						   Form("prfShiftQxCentEtaMinus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
      prfShiftQyCentEtaMinus[ii][jj] = new TProfile(Form("prfShiftQyCentEtaMinus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),
						    Form("prfShiftQyCentEtaMinus_harmonic_%i",ShiftConstants::mEPHarmonic[ii]*(jj+1)),9,0,9);
  
    }
  }
}
//-----------------------------------------------------------------------------
Int_t StShiftCorrection::Init()
{
  //cout<<"Lomnitz : pre enter"<<endl;
   if (mFileOut)
   {
     // cout<<"Lomnitz if passed"<<endl;
      mFileOut->cd();

      // event plane and Q vector
      float PI = TMath::Pi();
   }
   return kStOk;
}
//----------------------------------------------------------------------------- 
Int_t StShiftCorrection::Finish()
{
  cout<<"StShiftCorrection::Finish()"<<endl;
  mFileOut->cd();
  //
  for( int ii = 0 ; ii<ShiftConstants::mNoEP; ++ii){
    prfQxCentEtaPlus[ii]->Write();
    prfQyCentEtaPlus[ii]->Write();
    prfQxCentEtaMinus[ii]->Write();
    prfQxCentEtaMinus[ii]->Write();
    for(int jj = 0 ; jj < ShiftConstants::mShiftMaxHarmonic; ++jj){
      prfShiftQxCent[ii][jj]->Write();
      prfShiftQyCent[ii][jj]->Write();
      //
      prfShiftQxCentEtaPlus[ii][jj] ->Write();
      prfShiftQyCentEtaPlus[ii][jj] ->Write();
      prfShiftQxCentEtaMinus[ii][jj] ->Write();
      prfShiftQyCentEtaMinus[ii][jj] ->Write();
    }
  }
  //  mFileOut->Write();
  //  mFileOut->Close();

  return kStOK;
}
//-----------------------------------------------------------------------------
void StShiftCorrection::setFileOut(TFile* fileOut)
{
   mFileOut = fileOut;
}
//-----------------------------------------------------------------------------
Int_t StShiftCorrection::Make()
{
   if (!mPicoDstMaker)
   {
      LOG_ERROR << " No PicoDstMaker! Skip! " << endm;
      return kStErr;
   }

   mPicoDst = mPicoDstMaker->picoDst();
   if (!mPicoDst)
   {
      LOG_ERROR << " No PicoDst! Skip! " << endm;
      return kStErr;
   }

   mPicoEvent = (StPicoEvent*)mPicoDst->event();
   if (!mPicoEvent)
   {
      LOG_ERROR << "Error opening picoDst Event, skip!" << endm;
      return kStErr;
   }
   if (mRunNumber != mPicoEvent->runId()) getRunInfo(mPicoEvent->runId());
   else mAcceptQvectorFile = true;

   getEventInfo();//get event info

   if (mAcceptQvectorFile && mAcceptQvectorFiletmp)
   {
     bool bad_run = calculateEventPlane();
     if( !bad_run )
       calculateShiftCorrection();
   }

   return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StShiftCorrection::getEventInfo()
{

   //Remove bad vertices
   mVertexPos = mPicoEvent->primaryVertex();

   mgrefmultCorrUtil->init(mPicoDst->event()->runId());
   mgrefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mVertexPos.z(), mPicoDst->event()->ZDCx()) ;
   mCent  = mgrefmultCorrUtil->getCentralityBin9();

   mAcceptEvent = false;

   mBField = mPicoEvent->bField();

   bool isVPDMB5 = kFALSE;
   for (int i = 0; i < ShiftConstants::nTrig; i++)
   {
     if ( mPicoEvent->isTrigger(ShiftConstants::mTriggerId[i]) ) isVPDMB5 = kTRUE ;  //Select MB trigger
   }
   if (!(isVPDMB5))
   {
      return;
   }

   if (TMath::Abs(mVertexPos.z()) > ShiftConstants::vzMax) return;
   if (TMath::Abs(mVertexPos.z() - mPicoEvent->vzVpd()) > ShiftConstants::deltaVzMax) return;
   if (mCent < 0 || mCent > 9) return;

   mAcceptEvent = true;
}

void StShiftCorrection::getRunInfo(int const runNumber)
{
   mRunNumber = runNumber;

   char fileName[256];
   sprintf(fileName, "%s/%i.qVector.root", ShiftConstants::qVectorRunDir.Data(), mRunNumber);
   TFile* fQVector = new TFile(fileName);
   if(fQVector->IsZombie())
     {
       int dayNumber = mRunNumber%1000000/1000;
       sprintf(fileName, "%s/%03d.qVector.root", ShiftConstants::qVectorDayDir.Data(), dayNumber);
       delete fQVector;
       fQVector = new TFile(fileName);
       if(fQVector->IsZombie())
	 {
	   cout<<"can not load run or day qVector file: "<<mRunNumber<<endl;
	   return;
	 }
     }
   cout << "load qVector file: " << fileName << endl;

   for( int ii = 0; ii < ShiftConstants::mNoEP ; ++ii){
     int mHarmonic = ShiftConstants::mEPHarmonic[ii];
     fQVector->GetObject(Form("prfQxCentEtaPlus_v%i",mHarmonic), prfQxCentEtaPlus[ii]);
     if (!prfQxCentEtaPlus)
       {
	 LOG_INFO << "StShiftCorrection::THistograms and TProiles NOT found! shoudl check the files From HaoQiu" << endm;
	 mAcceptQvectorFile = false;
	 mAcceptQvectorFiletmp = false;
	 return;
       }
     else
       {
	 mAcceptQvectorFile = true;
	 mAcceptQvectorFiletmp = true;
       }
     
     prfQxCentEtaPlus[ii] = (TProfile*)fQVector->Get(Form("prfQxCentEtaPlus_v%i",mHarmonic))->Clone(Form("prfQxCentEtaPlus_v%i",mHarmonic));
     prfQyCentEtaPlus[ii] = (TProfile*)fQVector->Get(Form("prfQyCentEtaPlus_v%i",mHarmonic))->Clone(Form("prfQyCentEtaPlus_v%i",mHarmonic));
     prfQxCentEtaMinus[ii] = (TProfile*)fQVector->Get(Form("prfQxCentEtaMinus_v%i",mHarmonic))->Clone(Form("prfQxCentEtaMinus_v%i",mHarmonic));
     prfQyCentEtaMinus[ii] = (TProfile*)fQVector->Get(Form("prfQyCentEtaMinus_v%i",mHarmonic))->Clone(Form("prfQyCentEtaMinus_v%i",mHarmonic));

     prfQxCentEtaPlus[ii]->SetDirectory(0);
     prfQyCentEtaPlus[ii]->SetDirectory(0);
     prfQxCentEtaMinus[ii]->SetDirectory(0);
     prfQyCentEtaMinus[ii]->SetDirectory(0);
   }
   fQVector->Close();
   delete fQVector;
}
/*----------------------------------------------------------------------------------------------------------------------*/
bool StShiftCorrection::calculateEventPlane(){
   // track loop
   float Qx[ShiftConstants::mNoEP] = {0.};
   float Qy[ShiftConstants::mNoEP] = {0.};
   float QxEtaPlus[ShiftConstants::mNoEP] = {0.};
   float QyEtaPlus[ShiftConstants::mNoEP] = {0.};
   float QxEtaMinus[ShiftConstants::mNoEP] = {0.};
   float QyEtaMinus[ShiftConstants::mNoEP] = {0.};

   float vertexZ = mVertexPos.z();
   for (unsigned short iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }
      if (picoTrack->nHitsFit() < ShiftConstants::nHitsFitMin) continue;
      if(1.*picoTrack->nHitsFit()/picoTrack->nHitsMax() < ShiftConstants::mNHitsFitRatioMin) continue;

      StPhysicalHelix helix = picoTrack->dcaGeometry().helix();
      float dca = helix.geometricSignedDistance(mVertexPos);

      if (TMath::Abs(dca) > ShiftConstants::mDcaMax) continue;

      float pathLengthToPrimaryVertex = helix.pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix.momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      //SL16d
      double eta, pt, phi;
      if( picoTrack->pMom().perp() > 0 ){
	eta =picoTrack->pMom().pseudoRapidity();
	pt =picoTrack->pMom().perp();
	phi =picoTrack->pMom().phi();
      }
      else{
#ifdef only__primary
	continue;
#endif
	pt = momentum.perp();
	eta = momentum.pseudoRapidity();
	phi = momentum.phi();
      }

      if (fabs(eta) > ShiftConstants::mEtaMax) continue;
      if (pt < ShiftConstants::mPtMin || pt > ShiftConstants::mPtMax) continue;

      for(int iEP = 0 ; iEP < ShiftConstants::mNoEP; ++iEP){
	int const mHarmonic = ShiftConstants::mEPHarmonic[iEP];
	
	float qx = cos(mHarmonic * phi) * pt;
	float qy = sin(mHarmonic * phi) * pt;
	
	if (eta > 0)
	  {
	    qx -= prfQxCentEtaPlus[iEP]->GetBinContent(mCent + 1);
	    qy -= prfQyCentEtaPlus[iEP]->GetBinContent(mCent + 1);
	  }
	else
	  {
	    qx -= prfQxCentEtaMinus[iEP]->GetBinContent(mCent + 1);
	    qy -= prfQyCentEtaMinus[iEP]->GetBinContent(mCent + 1);
	  }
      
	Qx[iEP] += qx;
	Qy[iEP] += qy;
      
      if (eta > 0.05)
	{
	  QxEtaPlus[iEP] += qx;
	  QyEtaPlus[iEP] += qy;
	}
      if(eta < -0.05)
	{
	  QxEtaMinus[iEP] += qx;
	  QyEtaMinus[iEP] += qy;
	}
      
      //      iTrackForEventPlane++;
      }//loop thru picoTracks
   }
//assert(iTrackForEventPlane == nTracksForEventPlane);
   
   for(int iEP = 0 ; iEP < ShiftConstants::mNoEP; ++iEP){
     TVector2 mQ(Qx[iEP], Qy[iEP]);
     
     TVector2 mQEtaPlus(QxEtaPlus[iEP], QyEtaPlus[iEP]);
     TVector2 mQEtaMinus(QxEtaMinus[iEP], QyEtaMinus[iEP]);
     if (mQ.Mod2() == 0 ||  mQEtaPlus.Mod2() || mQEtaMinus.Mod2() == 0)
       {
	 return true;
       }
     
     mEventPlane[iEP] = mQ.Phi()/ShiftConstants::mEPHarmonic[iEP];
     mEventPlaneEtaPlus[iEP] = mQEtaPlus.Phi()/ShiftConstants::mEPHarmonic[iEP];
     mEventPlaneEtaMinus[iEP] = mQEtaMinus.Phi()/ShiftConstants::mEPHarmonic[iEP];
   }
   return false;
}
/*----------------------------------------------------------------------------------------------------------------------*/
void StShiftCorrection::calculateShiftCorrection()
{
  for(int iEP = 0 ; iEP < ShiftConstants::mNoEP; ++iEP){
    for(int iHarmonic = 0 ; iHarmonic < ShiftConstants::mShiftMaxHarmonic; ++iHarmonic){
      //      cout<<"Lomnitz"<<iEP<<" "<<iHarmonic<<endl;
      float const qx = cos(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlane[iEP]);
      float const qy = sin(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlane[iEP]);
      //
      float const qx_etaplus = cos(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlaneEtaPlus[iEP]);
      float const qy_etaplus = sin(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlaneEtaPlus[iEP]);
      float const qx_etaminus = cos(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlaneEtaMinus[iEP]);
      float const qy_etaminus = sin(ShiftConstants::mEPHarmonic[iEP] * (iHarmonic+1) * mEventPlaneEtaMinus[iEP]);

      prfShiftQxCent[iEP][iHarmonic]->Fill(mCent,qx);
      prfShiftQyCent[iEP][iHarmonic]->Fill(mCent,qy);
      //
      prfShiftQxCentEtaPlus[iEP][iHarmonic]->Fill(mCent,qx_etaplus);
      prfShiftQyCentEtaPlus[iEP][iHarmonic]->Fill(mCent,qy_etaplus);
      //
      prfShiftQxCentEtaMinus[iEP][iHarmonic]->Fill(mCent,qx_etaminus);
      prfShiftQyCentEtaMinus[iEP][iHarmonic]->Fill(mCent,qy_etaminus);
    }
  }
}
