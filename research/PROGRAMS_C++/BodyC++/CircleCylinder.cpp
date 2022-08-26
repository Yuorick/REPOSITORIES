//---------------------------------------------------------------------------


#pragma hdrstop

#include "CircleCylinder.h"
#include "Plane.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>




  TCircleCylinder ::TCircleCylinder()
{

	mR = 0. ;
	mH = 0.;
  //	mM = 0.;


}
TCircleCylinder ::TCircleCylinder(const TPlane Plane, const double R, const double H,const double M) :TSimpleBody_3D( Plane, M)
{
 mR = R;
 mH = H;
// mM = M;

}

// конструктор копирования
 TCircleCylinder ::TCircleCylinder (const TCircleCylinder &R)
  {
	mR  = R.mR;
	mH = R.mH;
	mM = R.mM   ;
	mPlane = R.mPlane;
 }
 // оператор присваивания
 TCircleCylinder TCircleCylinder::operator=(TCircleCylinder  R)
 {
	mR  = R.mR;
	mH = R.mH;
	mM = R.mM   ;
	mPlane = R.mPlane;
	return *this ;
  }


// вычисление объема
double TCircleCylinder::calcCapacity()
{
 return  M_PI * mH * mR * mR ;
}

// вычисление координат центра тяжести
void TCircleCylinder::calcCentreOfGravity(double *arrCentreGrav)
{
  arrCentreGrav [0] = mH /2.;
  arrCentreGrav [1] = 0.;
  arrCentreGrav [2] = 0.;
}

// вычисление тензора моментов инерции
// Динник Справочник по технической механике, стр 236
// это диагональная матрица

void TCircleCylinder::calcInertiaMtrx(double *arrInertMtrx)
{
  // 1. высисление момента инерции относительно оси Y
  memset(arrInertMtrx, 0, 9 * sizeof(double));
  arrInertMtrx[0] =  mR*mR / 2.;
  ////
  arrInertMtrx[4] = (3. * mR * mR + mH * mH)/ 12.;
  arrInertMtrx[8] = arrInertMtrx[4];
}


#pragma package(smart_init)
