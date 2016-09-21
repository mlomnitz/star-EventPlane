/* **************************************************
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu         (hqiu@lbl.gov)
 *            **Guannan Xie     (guannanxie@lbl.gov)
 *
 *            ** code maintainer
 *
 * **************************************************
 */


#ifndef STAR_StShiftCorrection
#define STAR_StShiftCorrection

#include "StMaker.h"
#include "TVector2.h"
#include "StThreeVectorF.hh"

#include "StShiftConstants.h"

class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StRefMultCorr;
class TH1I;
class TH1F;
class TH2F;
class TH3F;
class THn;
class TProfile;

const int maxNTracks = 20000;

class StShiftCorrection : public StMaker
{
public:
  StShiftCorrection(const char* name, StPicoDstMaker* picoMaker, StRefMultCorr* grefmultCorrUtil);

   Int_t Init();
   virtual Int_t Make();
   Int_t Finish();
   void setFileOut(TFile* fileOut);

   int   getRunId() const;
   bool getAcceptEvent() const;

private:
   void getEventInfo();
   void getRunInfo(int runNumber);
   void calculateShiftCorrection();
   bool calculateEventPlane();

   StPicoDstMaker* mPicoDstMaker;
   StPicoDst*      mPicoDst;
   StPicoEvent*    mPicoEvent;
   StRefMultCorr* mgrefmultCorrUtil;

   bool   mAcceptEvent;
   bool   mAcceptQvectorFile;
   bool   mAcceptQvectorFiletmp;

   int         mCent;
   int         mRunNumber;
   float       mBField;
   StThreeVectorF mVertexPos;

   float mEventPlane[ShiftConstants::mNoEP];
   float mEventPlaneEtaPlus[ShiftConstants::mNoEP];
   float mEventPlaneEtaMinus[ShiftConstants::mNoEP];

   TFile* mFileOut;

   TProfile* prfQxCentEtaPlus[ShiftConstants::mNoEP];
   TProfile* prfQyCentEtaPlus[ShiftConstants::mNoEP];
   TProfile* prfQxCentEtaMinus[ShiftConstants::mNoEP];
   TProfile* prfQyCentEtaMinus[ShiftConstants::mNoEP];
   //
   TProfile* prfShiftQxCent[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];
   TProfile* prfShiftQyCent[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];
   //
   TProfile* prfShiftQxCentEtaPlus[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];
   TProfile* prfShiftQyCentEtaPlus[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];
   TProfile* prfShiftQxCentEtaMinus[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];
   TProfile* prfShiftQyCentEtaMinus[ShiftConstants::mNoEP][ShiftConstants::mShiftMaxHarmonic];

   ClassDef(StShiftCorrection, 0)
};

#endif
