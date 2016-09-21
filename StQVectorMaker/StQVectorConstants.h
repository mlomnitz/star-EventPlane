#ifndef __qVectorConst__
#define __qVectorConst__
#include "StQVectorConstants.h"

namespace qVectorConst{
      //Event Cuts 
  float const mVzMax = 6.0;
  float const mDeltaVzMax = 3.0;
  //Track Cuts
  float const mNHitsFitMin = 15;
  float const mNHitsFitRatioMin = 0.52;
  float const mEtaMax = 1.0;
  float const mPtMin = 0.2;
  float const mPtMax = 2.;
  float const mDcaMax = 1.0;
  //
  int const mNoEP = 2;
  int const mEPHarmonic[mNoEP] = {2,3};
}
#endif
