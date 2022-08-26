//---------------------------------------------------------------------------


#pragma hdrstop

#include "InitTargData.h"
#include <math.h>


TInitTargData::TInitTargData()
{
	//-  угол пеленга цели в ГСК
	mBearing = M_PI/2.;
	// угол курса цели в ГСК
	mTargCourse = 1.5 * M_PI;
	// угол между ветором скорости цели и напрвлением в зенит
	mTargZenitAng = M_PI/2. ;
	//- скорость
	mV  = 1000. ;
	// дальность,
	mR =20000. ;
	// высота
	mH = 1000. ;
	// время
	mT = 0. ;

}

//---------------------------------------------------------------------------
  __fastcall TInitTargData::~TInitTargData()
{

}

// конструктор копирования
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
 // оператор присваивания
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

  // парам конструктор1
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
