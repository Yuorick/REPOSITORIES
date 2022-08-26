//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "Sector.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 TSector ::TSector()
{
	mPntCentre = TURPointXY(0., 0.);
	mR = 3. ;
	mFi0 = M_PI ;
	mFi1 =  M_PI/ 2.;
}

 // парам констр
TSector :: TSector( const TURPointXY PntCentre, const double R,  const double Fi0,  const double Fi1 )

 {
   mPntCentre = PntCentre ;
   mR = R ;
   mFi0 = Fi0 ;
   mFi1 = Fi1;
 }





 // оператор присваивания
 TSector &TSector::operator=(const TSector  &R)
 {
	mPntCentre = R.mPntCentre ;
	mR = R.mR ;
	mFi0 = R.mFi0 ;
	mFi1 = R.mFi1 ;
	return *this ;
 }

 // конструктор копирования
 TSector::TSector (const TSector &R)
 {
	mPntCentre = R.mPntCentre ;
	mR = R.mR ;
	mFi0 = R.mFi0 ;
	mFi1 = R.mFi1 ;
 }

void TSector::ShowMe(wchar_t *FileName)
{
   int iN = 1200;
   TURPolyLine line(*this,iN);
   line.WriteSetSHPFiles(FileName,&line, 1 ) ;
}



#pragma package(smart_init)
