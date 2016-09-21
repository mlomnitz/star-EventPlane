#define only__primary
#include "StQVectorMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "PhysicalConstants.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TVector2.h"

#include "StQVectorHists.h"
ClassImp(StQVectorMaker)

using namespace qVectorConst;
//-----------------------------------------------------------------------------
StQVectorMaker::StQVectorMaker(const char* name, StPicoDstMaker *picoMaker)
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
}

//-----------------------------------------------------------------------------
StQVectorMaker::~StQVectorMaker()
{ /*  */ }

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Init()
{

    mAcceptEvent = false;

    mRefMultCorr = new StRefMultCorr("grefmult");

    mFileOut = new TFile(mOutputName, "recreate");

    // event level QA
    hVzVpdVz = new TH2F("hVzVpdVz","hVzVpdVz",200,-100,100,200,-100,100);
    hVzDiff  = new TH1F("hVzDiff","hVzDiff",500,-100,100);
    hVxy = new TH2F("hVxy","hVxy",500,-1,1,500,-1,1);
    hRefMult = new TH1I("hRefMult","hRefMult",1000,0,1000);
    hGRefMult = new TH1I("hGRefMult","hGRefMult",1000,0,1000);
    hTrigger = new TH1I("hTrigger","hTrigger",32,0,32);
    hCentrality = new TH1I("hCentrality", "hCentrality", 9, 0, 9);

    // track level QA
    hNHitsFit = new TH1I("hNHitsFit", "hNHitsFit", 50, 0, 50);
    hDca = new TH1F("hDca", "hDca", 20000, -10, 10);
    hEta = new TH1F("hEta", "hEta", 30, -1.5, 1.5);
    hPt = new TH1F("hPt", "hPt", 200, 0, 10);
    
    // Phi dist
    float PI = TMath::Pi();

    hPhiCentEtaPlusZPlus = new TH2F("hPhiCentEtaPlusZPlus","hPhiCentEtaPlusZPlus",9,0,9,120,-PI,PI);
    hPhiCentEtaPlusZMinus = new TH2F("hPhiCentEtaPlusZMinus","hPhiCentEtaPlusZMinus",9,0,9,120,-PI,PI);
    hPhiCentEtaMinusZPlus = new TH2F("hPhiCentEtaMinusZPlus","hPhiCentEtaMinusZPlus",9,0,9,120,-PI,PI);
    hPhiCentEtaMinusZMinus = new TH2F("hPhiCentEtaMinusZMinus","hPhiCentEtaMinusZMinus",9,0,9,120,-PI,PI);

    for(int ii = 0; ii<mNoEP; ++ii){
      eventPlane[ii] = new StQVectorHists(mEPHarmonic[ii]);
    }

    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Finish() {
    mFileOut->Write();
    return kStOK;
}

//-----------------------------------------------------------------------------
void StQVectorMaker::Clear(Option_t *opt) {
}

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Make() {
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    getEventInfo();//get event info
    if(mAcceptEvent){
        getTrackInfo();
    }

    return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::setOutputName(Char_t* dir, Char_t* file)
{
    TString dirName(dir);
    TString fileName(file);
    mOutputName = dirName+"/"+fileName+".qVector.root";	
    
    return;
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::getEventInfo()
{
    mAcceptEvent = false;
    
    if(!mPicoDst) return;
    
    //Load event
    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent){
        cerr<<"Error opening picoDst Event, skip!"<<endl;
        return;
    }
    
    for(int i=0; i<32; i++)
      if( mPicoEvent->isTrigger(i) )
            hTrigger->Fill(i);

    if( !isMinBiasTrigger() )
      return;
    
    //Remove bad vertices
    mVertexPos = mPicoEvent->primaryVertex();
    
    hVzVpdVz->Fill(mVertexPos.z(), mPicoEvent->vzVpd());
    hVzDiff->Fill(mPicoEvent->vzVpd() - mVertexPos.z());
    hVxy->Fill(mVertexPos.y(), mVertexPos.x());

    if(TMath::Abs(mVertexPos.z()) > mVzMax) return;
    if(TMath::Abs(mVertexPos.z() - mPicoEvent->vzVpd()) > mDeltaVzMax) return;

    hRefMult->Fill(mPicoEvent->refMult());
    hGRefMult->Fill(mPicoEvent->grefMult());

    if(!mRefMultCorr) {  
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return;
    } 
    mRefMultCorr->init(mPicoDst->event()->runId());
    mRefMultCorr->initEvent(mPicoDst->event()->grefMult(),mVertexPos.z(),mPicoDst->event()->ZDCx()) ;
    mCent  = mRefMultCorr->getCentralityBin9();

    if (mCent < 0 || mCent > 8) return;
    
    hCentrality->Fill(mCent);

    mAcceptEvent = true;

    mBField = mPicoEvent->bField();
}

/*----------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::getTrackInfo()
{
    float vertexZ = mPicoEvent->primaryVertex().z();
    //float Qx=0., Qy=0.;
    float Qx[mNoEP] = {0};
    float Qy[mNoEP] = {0};

    //Load tracks for default vertex index
    for(int iTrack=0; iTrack<mPicoDst->numberOfTracks(); iTrack++)
        {
            StPicoTrack* picoTrack =(StPicoTrack*) mPicoDst->track(iTrack);
            if(!picoTrack){
                break;
            }
            
            hNHitsFit->Fill(picoTrack->nHitsFit());
            if(picoTrack->nHitsFit() <= mNHitsFitMin) continue;
	    if(1.*picoTrack->nHitsFit()/picoTrack->nHitsMax() < mNHitsFitRatioMin) continue;

	    StPhysicalHelix helix = picoTrack->dcaGeometry().helix();
	    float dca = helix.geometricSignedDistance(mVertexPos);
            hDca->Fill(dca);
            if(TMath::Abs(dca) > mDcaMax) continue;
            
	    float pathLengthToPrimaryVertex =helix.pathLength(mVertexPos.x(), mVertexPos.y());
	    StThreeVectorF momentum = helix.momentumAt(pathLengthToPrimaryVertex, mBField*kilogauss);
	    double pt, eta, phi;
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
	    if( phi < 0 ) phi+=2.*TMath::Pi();
            hEta->Fill(eta);
            hPt->Fill(pt);
            if(fabs(eta) > mEtaMax) continue;
            if(pt<mPtMin || pt>mPtMax) continue;

            if(eta>0 && vertexZ>0) hPhiCentEtaPlusZPlus->Fill(mCent, phi);
            if(eta>0 && vertexZ<0) hPhiCentEtaPlusZMinus->Fill(mCent, phi);
            if(eta<0 && vertexZ>0) hPhiCentEtaMinusZPlus->Fill(mCent, phi);
            if(eta<0 && vertexZ<0) hPhiCentEtaMinusZMinus->Fill(mCent, phi);
	    
	    for(int iEP = 0; iEP < mNoEP; ++iEP){
	      float qx = cos(mEPHarmonic[iEP]*phi)*pt;
	      float qy = sin(mEPHarmonic[iEP]*phi)*pt;
	     
	      eventPlane[iEP]->addTrack(mCent, eta, qx, qy);
	      
	      Qx[iEP] += qx;
	      Qy[iEP] += qy;
	    }
        }//loop thru picoTracks
    for(int iEP = 0 ; iEP  < mNoEP; ++iEP)
      eventPlane[iEP]->addEventPlane(mCent, Qx[iEP], Qy[iEP]);
}
bool StQVectorMaker::isMinBiasTrigger() const 
{
  int mTriggerId[5] = {450050, 450060,
				  450005, 450015,
				  450025 };
  for(int ii = 0; ii<5; ++ii){
    if( mPicoEvent->isTrigger(mTriggerId[ii]) )
      return true;
  }
  return false;
}
