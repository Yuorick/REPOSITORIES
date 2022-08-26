//---------------------------------------------------------------------------


#pragma hdrstop

#include "SimpleDifDgrm.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
#include "MatrixProccess.h"

_fastcall  TSimpleDifDgrm ::TSimpleDifDgrm()
{
  mTang = 60.;;
 // угол сканирования разностной диаграммы
  mScnDif = -5. * 3000./ M_PI;

}
// Конструктор копирования
__fastcall  TSimpleDifDgrm::TSimpleDifDgrm (const TSimpleDifDgrm &R)
 {
	mTang =  R.mTang;
	// угол сканирования разностной диаграммы
	mScnDif = R.mScnDif;
	// понижающий коэфф суммарнй диаграммы

 }
// оператор присваивания
  TSimpleDifDgrm TSimpleDifDgrm::operator=(TSimpleDifDgrm  R)
 {
	mTang =  R.mTang;
	// угол сканирования разностной диаграммы
	mScnDif = R.mScnDif;

	 return *this ;
 }

 // парам констр
 __fastcall TSimpleDifDgrm::TSimpleDifDgrm(const double Tang, const double ScnDif )
 {
	 mTang = Tang ;
	 mScnDif = ScnDif;

 }


#pragma package(smart_init)
