//---------------------------------------------------------------------------


#pragma hdrstop

#include "RectWaveGuide.h"
#include <math.h>


//--------------------------------------------------------------------------------------
 TRectWaveGuide ::TRectWaveGuide()
{
	ma= 1.35;
	mb =   4. ;

}

 // парам констр
TRectWaveGuide :: TRectWaveGuide(  const double a,  const double b )

 {
   ma = a;
   mb = b;
 }





 // оператор присваивания
 TRectWaveGuide TRectWaveGuide::operator=(TRectWaveGuide  R)
 {
	ma = R.ma ;
	mb = R.mb ;
	return *this ;
 }

 // конструктор копирования
 TRectWaveGuide::TRectWaveGuide (const TRectWaveGuide &R)
 {
	ma = R.ma ;
	mb = R.mb ;
 }
/*
void TRectWaveGuide::ShowMe(wchar_t *FileName)
{
   int iN = 1200;
   TURPolyLine line(*this,iN);
   line.WriteSetSHPFiles(FileName,&line, 1 ) ;
}
 */
double TRectWaveGuide::fncDiagr(const double VAlLambda, const double VAlTetta)
{
  if (fabs(VAlTetta) < 0.0000000001)
  {
   return 1.;
  }
  double temp = M_PI * mb / VAlLambda;
  return exp(- temp * sin(VAlTetta/ 2.)* sin(VAlTetta/ 2.))* sqrt(sin(temp * sin(VAlTetta)) / (temp * sin(VAlTetta)));
}

#pragma package(smart_init)
