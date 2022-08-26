//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "Environment.h"

//---------------------------------------------------------------------------

__fastcall TEnvironment::TEnvironment()
{
	// �������� �� ���� �� ����� Significance Wave Height (SWH). �� 0 �� 9 ������
	 mBallWave = 0.;
//  �������� �����  ��������������  �/�
	 mWind_V = 0. ;

//  �������� �����  ������������  �/�
	 mWind_VertV = 0. ;
	 // ����������� ������ ���� ������
	 mWind_Alf = 0. ;

}

//---------------------------------------------------------------------------


// ����������� �����������
  __fastcall TEnvironment ::TEnvironment (const TEnvironment &R)
 {
   mBallWave = R.mBallWave;
   mWind_V = R.mWind_V ;
   mWind_Alf = R.mWind_Alf;
   mWind_VertV = R.mWind_VertV;

 }
 // �������� ������������
 TEnvironment &TEnvironment::operator=(const TEnvironment  &R)
 {
   mBallWave = R.mBallWave;
   mWind_V = R.mWind_V ;
   mWind_Alf = R.mWind_Alf;
   mWind_VertV = R.mWind_VertV;
	return *this ;
 }


  // ����� �����������
  __fastcall TEnvironment::TEnvironment (const int BallWave, const double Wind_V)
 {
	mBallWave  = BallWave ;
	mWind_V  = Wind_V  ;
	mWind_VertV = 0.;

 }

  // ����� �����������
  __fastcall TEnvironment::TEnvironment (const double Wind_V, const double Wind_Alf, const double  Wind_VertV)
 {
	mWind_V  = Wind_V  ;
	mWind_Alf = Wind_Alf;
	mBallWave  = int (Wind_V / 20. * 9. +0.0001)  ;
	mWind_VertV = Wind_VertV;

 }
#pragma package(smart_init)
