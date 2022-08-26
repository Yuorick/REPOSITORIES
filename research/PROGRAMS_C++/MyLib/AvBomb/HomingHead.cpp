//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <string.h>
#include "HomingHead.h"
THomingHead ::THomingHead()
{
	// ����
 mFi = 14. * M_PI / 180. ;
 // ���������
 mR= 5000.;
  // ����� �� ���������, �
 mSigFl_R = 3./2.;
 // ����� �� ����, ����
 mSigFl_U = 0.25 / 180. * M_PI/ 3.;
 m_h = 0.0002;
}



  // ����� ������  1
THomingHead :: THomingHead( const long  double Fi, const long  double R
   ,  const long  double SigFl_R, long  double SigFl_U , long  double h   )

 {
   mFi = Fi ;
   mR =R;
   mSigFl_R = SigFl_R;
   mSigFl_U =SigFl_U ;
   m_h = h ;
 }





 // �������� ������������
 THomingHead THomingHead::operator=(THomingHead  R)
 {
   mFi = R.mFi ;
   mR = R.mR;
   mFi = R.mFi ;
   mR = R.mR;
   mSigFl_R = R.mSigFl_R;
   mSigFl_U = R.mSigFl_U ;
   m_h = R.m_h ;
   return *this ;
 }

 // ����������� �����������
 THomingHead::THomingHead (const THomingHead &R)
 {
   mFi = R.mFi ;
   mR = R.mR;
   mFi = R.mFi ;
   mR = R.mR;
   mSigFl_R = R.mSigFl_R;
   mSigFl_U = R.mSigFl_U ;
   m_h = R.m_h ;
 }


#pragma package(smart_init)
