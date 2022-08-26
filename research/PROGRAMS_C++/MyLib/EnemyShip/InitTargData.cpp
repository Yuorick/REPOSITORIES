//---------------------------------------------------------------------------


#pragma hdrstop

#include "InitTargData.h"
#include <math.h>


TInitTargData::TInitTargData()
{
	//-  ���� ������� ���� � ���
	mBearing = M_PI/2.;
	// ���� ����� ���� � ���
	mTargCourse = 1.5 * M_PI;
	// ���� ����� ������� �������� ���� � ����������� � �����
	mTargZenitAng = M_PI/2. ;
	//- ��������
	mV  = 1000. ;
	// ���������,
	mR =20000. ;
	// ������
	mH = 1000. ;
	// �����
	mT = 0. ;

}

//---------------------------------------------------------------------------
  __fastcall TInitTargData::~TInitTargData()
{

}

// ����������� �����������
 TInitTargData ::TInitTargData (const TInitTargData &R)
 {
	mBearing  =  R.mBearing ;
	mTargCourse = R.mTargCourse;
	mTargZenitAng = R.mTargZenitAng ;
	mV   = R.mV ;
	mR = R.mR ;
	mH = R.mH ;
	mT = R.mT ;

 }
 // �������� ������������
 TInitTargData TInitTargData::operator=(TInitTargData  R)
 {
	mBearing  =  R.mBearing ;
	mTargCourse = R.mTargCourse;
	mTargZenitAng = R.mTargZenitAng ;
	mV   = R.mV ;
	mR = R.mR ;
	mH = R.mH ;
	mT = R.mT ;



	return *this ;
 }

  // ����� �����������1
 TInitTargData::TInitTargData (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
	 ,const double  R,const double H, const double T)
 {
   mBearing  = Bearing ;
	mTargCourse = TargCourse ;
	mTargZenitAng = TargZenitAng;
	mV   = V ;
	mR = R ;
	mH = H ;
	mT = T ;

 }



//---------------------------------------------------------------------------

#pragma package(smart_init)
