//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
//#include <stdlib.h>
#include "Zamer.h"

__fastcall  TZamer::TZamer()
{
 memset(marrMeas, 0, 3 * sizeof(double));
 memset(marrCorr, 0, 9 * sizeof(double));
 mT = 0.;

}
// Конструктор копирования
__fastcall  TZamer::TZamer (const TZamer &R2)
 {
	 memcpy(marrMeas, R2.marrMeas , 3 * sizeof(double));
	 memcpy(marrCorr, R2.marrCorr , 9 * sizeof(double));
	 mT = R2.mT;
 }
 // оператор присваивания
	TZamer &TZamer::operator=(const TZamer  &R2)
 {
	 memcpy(marrMeas, R2.marrMeas , 3 * sizeof(double));
	 memcpy(marrCorr, R2.marrCorr , 9 * sizeof(double));
	 mT = R2.mT;

	 return *this ;
 }

 // парам констр 1
 __fastcall TZamer::TZamer ( double *arrMeas,double *arrCorr, const double VAlT)
 {
	 memcpy(marrMeas, arrMeas , 3 * sizeof(double));
	 memcpy(marrCorr, arrCorr , 9 * sizeof(double));
	 mT = VAlT;
 }
#pragma package(smart_init)
