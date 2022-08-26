//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "Ebomb.h"
#include <string.h>
#include "UrPointXY.h"
#include "URPolyLine.h"


TEbomb ::TEbomb()
{
	//  угол раствора диаграммы
 mFi = 3. * M_PI / 180. ;
 // дальность диаграммы
 mR= 1000.;
 // требуемый поток мощности
 mW = 140.;

}

 // парам констр
TEbomb :: TEbomb(const long  double Fi,const long  double R,const long  double W )

 {
   mFi = Fi ;
   mR =R;
   mW = W ;
 }





 // оператор присваивания
 TEbomb TEbomb::operator=(TEbomb  R)
 {
   mFi = R.mFi ;
   mR = R.mR;
   mW = R.mW ;

	return *this ;
 }

 // конструктор копирования
 TEbomb::TEbomb (const TEbomb &R)
 {
   mFi = R.mFi ;
   mR = R.mR;
   mW = R.mW ;
 }

 // расчет мощности

long double TEbomb::fncCalcPower ()
{
  long double valTelAng =  2. * M_PI * (1. - cosl(mFi / 2.));
  return mW * mR *mR * valTelAng *10000.;
}

void TEbomb::createGraphR_of_Fi (wchar_t *pwchOutFile, double valW, double valP)
{
  double valFi0 = 3.7 * M_PI/ 180.;
  double valFi1 = M_PI/ 2.;
  int N = 100;
  double step = (valFi1 - valFi0)/((double) N);
  TURPointXY *pntarr = new TURPointXY[N] ;
  for (int i = 0; i < N; i++)
  {
	double valFi = valFi0 + ((double)i)  *  step;
	pntarr[i].X =  valFi * 100.;
	pntarr[i].Y =  sqrt(valW* 100000./ valP/ 2./M_PI/ ( 1. - cos(valFi/2.)));

  }
	TURPolyLine  PolyLine( pntarr,N) ;
	delete [] pntarr;
	wchar_t wchFileName[300]  = {0};
	wcscpy(wchFileName,pwchOutFile );
	wcscat(wchFileName, L"\\GraphRDefeat_of_Fi.shp");

	PolyLine.WriteSetSHPFiles(wchFileName,&PolyLine, 1) ;
}

#pragma package(smart_init)
