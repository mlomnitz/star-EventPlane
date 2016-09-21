/* **************************************************
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu         (hqiu@lbl.gov)
 *            **Guannan Xie     (guannanxie@lbl.gov)
 *
 *            ** code maintainer
 *
 * **************************************************
 */


#ifndef StShiftConstants_H
#define StShiftConstants_H

#include "TString.h"

namespace ShiftConstants
{
  int const nTrig = 5;
  int const mTriggerId[nTrig] = {450050, 450060,
				 450005, 450015,
				 450025 };

  //SL16d prduction
  TString const qVectorRunDir ="/global/homes/m/mlomnitz/mlomnitz_projectdir/Run14_Reprod/recenter_v2_2/qVectorRun_vn_primary_only";
  TString const qVectorDayDir = "/global/homes/m/mlomnitz/mlomnitz_projectdir/Run14_Reprod/recenter_v2_2/qVectorDay_vn_primary_only";
   //Event Cuts
   float const vzMax = 6.0;
   float const deltaVzMax = 3.0;

   //Track Cuts
   int const nHitsFitMin = 15;
   float const mNHitsFitRatioMin = 0.52;
   //Track cuts for event plane
   float const mEtaMax = 1.0;
   float const mPtMin  = 0.2;
   float const mPtMax  = 2.0;
   float const mDcaMax = 1.0;
   //
   int const mNoEP = 2;
   int const mEPHarmonic[mNoEP] = {2,3};
   //
   int const mShiftMaxHarmonic = 2;
   
}
#endif
