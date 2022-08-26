//---------------------------------------------------------------------------


#pragma hdrstop
#include "Far.h"

 #include "Faceta.h"
//---------------------------------------------------------------------------

#include <vcl.h>


#include <math.h>
#include "Comp.h"
#include <stdio.h>
#include <stdlib.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

TFar::TFar()
{
 // количество модулей (фасет)
  m_N = 4;
 // расстояние между модулями
   m_D = 32.8;
 // длина волны
   mLambda = 3.15;;
 // массив фасет
  mFaceta = TFaceta();

}
// Конструктор копирования
TFar::TFar (const TFar &R2)
 {
 // количество модулей (фасет)
  m_N = R2.m_N;
 // расстояние между модулями
  m_D = R2.m_D;
 // длина волны
  mLambda = R2.mLambda;

  mFaceta = R2.mFaceta;

 }

// парам констр
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta)
 {
	 m_N = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;

 }

 // парам констр
 __fastcall TFar::TFar(const int N)
 {

 // количество модулей (фасет)
  m_N = N;
 // расстояние между модулями
   m_D = 32.8;
 // длина волны
   mLambda = 3.15;;
 // массив фасет
  mFaceta = TFaceta();

 }
// {

// }

 // оператор присваивания
 TFar TFar::operator=(TFar  R2)
{
 // количество модулей (фасет)
  m_N = R2.m_N;
 // расстояние между модулями
  m_D = R2.m_D;
 // длина волны
  mLambda = R2.mLambda;

  mFaceta = R2.mFaceta;

  return *this ;
}



#pragma package(smart_init)
