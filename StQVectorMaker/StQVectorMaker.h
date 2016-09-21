#ifndef STAR_StQVectorMaker
#define STAR_StQVectorMaker
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "StQVectorConstants.h"
class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StRefMultCorr;
class TH1I;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class StQVectorHists;

class StQVectorMaker : public StMaker {
  public:
    StQVectorMaker(const char *name, StPicoDstMaker *picoMaker);
    virtual ~StQVectorMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    void getEventInfo();
    void getTrackInfo();
    void setOutputName(Char_t* dir=".", Char_t* name="test");
    bool isMinBiasTrigger() const;
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent	   *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;

    bool	mAcceptEvent;
    int         mCent;
    float       mBField;
    StThreeVectorF mVertexPos;

    TString mOutputName;
    TFile* mFileOut;
        
    //event level qa
    TH2F*      hVzVpdVz;
    TH1F*      hVzDiff;
    TH2F*      hVxy;
    TH1I*      hRefMult;
    TH1I*      hGRefMult;
    TH1I*      hTrigger;
    TH1I*      hCentrality;

    //track level qa
    TH1I*      hNHitsFit;
    TH1F*      hDca;
    TH1F*      hEta;
    TH1F*      hPt;
    // Event Planes
    TH2F*      hPhiCentEtaPlusZPlus;
    TH2F*      hPhiCentEtaPlusZMinus;
    TH2F*      hPhiCentEtaMinusZPlus;
    TH2F*      hPhiCentEtaMinusZMinus;
    StQVectorHists* eventPlane[qVectorConst::mNoEP];


    ClassDef(StQVectorMaker, 1)
};

#endif
