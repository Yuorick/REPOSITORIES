//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "TargBearing0.h"


TTargBearing0::TTargBearing0()
{
	//-  угол пеленга цели в √— 
	mBearing = M_PI/2.;
	// угол курса цели в √— 
	mTargCourse = 1.5 * M_PI;
	// угол между ветором скорости цели и напрвлением в зенит
	mTargHorizAng = M_PI/2. ;
	//- скорость
	mV  = 1000. ;
	// дальность,
	mR =20000. ;
	// высота
	mH = 1000. ;
	// врем€
  //	mMy  = 1. ;
	//
   //	mTargType = PLANE;
	//
  //	mTargEPR = 1.;

}

//---------------------------------------------------------------------------
  __fastcall TTargBearing0::~TTargBearing0()
{

}

// конструктор копировани€
 TTargBearing0 ::TTargBearing0 (const TTargBearing0 &R)
 {
	mBearing  =  R.mBearing ;
	mTargCourse = R.mTargCourse;
	mTargHorizAng = R.mTargHorizAng ;
	mV   = R.mV ;
	mR = R.mR ;
	mH = R.mH ;
 //	mMy = R.mMy ;
  //	mTargType = R.mTargType;
   //	mTargEPR = R.mTargEPR;

 }
 // оператор присваивани€
 TTargBearing0 TTargBearing0::operator=(TTargBearing0  R)
 {
	mBearing  =  R.mBearing ;
	mTargCourse = R.mTargCourse;
	mTargHorizAng = R.mTargHorizAng ;
	mV   = R.mV ;
	mR = R.mR ;
	mH = R.mH ;
   //	mMy  = R.mMy  ;
  //	mTargType = R.mTargType;
  //	mTargEPR = R.mTargEPR;

	return *this ;
 }

  // парам конструктор1
 TTargBearing0::TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
	 ,const double  R,const double H)//, const double My), enumTargetType TargType )
 {
   mBearing  = Bearing ;
	mTargCourse = TargCourse ;
	mTargHorizAng = TargZenitAng;
	mV   = V ;
	mR = R ;
	mH = H ;
   //	mMy  = My  ;
  //	mTargType = TargType;

 }
 /*
 // парам конструктор1
 TTargBearing0::TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
	 ,const double  R,const double H, const double My, const double TargEPR , enumTargetType TargType )
 {
   mBearing  = Bearing ;
	mTargCourse = TargCourse ;
	mTargHorizAng = TargZenitAng;
	mV   = V ;
	mR = R ;
	mH = H ;
  //	mMy  = My  ;
  //	mTargType = TargType;
   //	mTargEPR = TargEPR;
 }
 */

 // расчет начальных координат цели
 // pX[6] - вектор координат цели (3 положени€, 3 скорости)
void TTargBearing0::raschet_nach_coord (double *pX)
{
  const double VAlE0=asin(mH / mR);
  // ¬ычисление начальных координат цели
  pX[0] = mR *cos(VAlE0)*sin(mBearing);
  pX[1] = (mR) *cos(VAlE0)*cos(mBearing);
 // if(mTargType == ABOVEWATER) pX[2] =5.;
 // else
  pX[2] = mR *sin(VAlE0);

  pX[3] =  mV * cos(mTargHorizAng) * sin (mTargCourse);
  pX[4] =  mV * cos(mTargHorizAng) * cos (mTargCourse);
  pX[5] =  mV * sin(mTargHorizAng) ;
}


#pragma package(smart_init)
