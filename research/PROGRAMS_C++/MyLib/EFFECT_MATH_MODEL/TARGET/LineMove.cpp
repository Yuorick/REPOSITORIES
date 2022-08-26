//---------------------------------------------------------------------------


#pragma hdrstop

#include "LineMove.h"

#include "LineMove.h"
#include "Traject.h"
#include <math.h>
#include  <string.h>

TLineMove::TLineMove()
{

	mT0 = 0. ;
	memset(marrVS, 0, 9 * sizeof(double)) ;

}

//---------------------------------------------------------------------------
  __fastcall TLineMove::~TLineMove()
{

}

// конструктор копирования
 TLineMove ::TLineMove (const TLineMove &R)
 {
  mT0  = R.mT0;
  memcpy( marrVS, R.marrVS, 9 * sizeof(double)) ;

 }
 // оператор присваивания
 TLineMove TLineMove::operator=(TLineMove  R)
 {
  mT0  = R.mT0;
  memcpy( marrVS, R.marrVS, 9 * sizeof(double)) ;
  return *this ;
 }

  // парам конструктор1
 TLineMove::TLineMove (const double T,  double *arrVS)
 {
	mT0  = T ;
	memcpy( marrVS, arrVS, 9 * sizeof(double)) ;

 }
// пролонгация вектора состояния на момент tExtr
void TLineMove::ExtrapolateVS (const double tExtr, double *arrVSTargExtr_GSK)
{
  for (int i = 0; i < 3 ;i++)
  {
	arrVSTargExtr_GSK[i] = marrVS[i] + (tExtr - mT0) * marrVS[3 + i] ;
	arrVSTargExtr_GSK [3 + i] =  marrVS[3 + i] ;
  }

}




#pragma package(smart_init)
