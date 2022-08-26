﻿//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "MatrixProccess.h"
#include "SimpleBody_3D.h"

__fastcall TSimpleBody_3D::TSimpleBody_3D()
{
	/*RecNumber = 0 ;
	RecLength = 0 ;
	Type = ShapeType::NullShape ;
	FigureColor = (TColor)0 ;*/
	 mPlane = TPlane();
	 mM = 0.; // масса
}
// конструктор копирования
 TSimpleBody_3D ::TSimpleBody_3D  (const TSimpleBody_3D &R)
  {
	mPlane  = R.mPlane;
	mM = 0.;
 }

  // оператор присваивания
  TSimpleBody_3D TSimpleBody_3D::operator=(TSimpleBody_3D  R)
  {
	mPlane  = R.mPlane;
	mM = 0.;

	return *this;
 }
//---------------------------------------------------------------------------
__fastcall TSimpleBody_3D::TSimpleBody_3D(const TPlane Plane, const double M)
{
	 mPlane = Plane;
	 mM = M;
}
//---------------------------------------------------------------------------

//__fastcall TSimpleBody_3D::~TSimpleBody_3D()
//{
//}
//---------------------------------------------------------------


// вычисление объема
double TSimpleBody_3D::calcCapacity()
{

}

// вычисление координат центра тяжести  в системе собственных координат простого тела
void TSimpleBody_3D::calcCentreOfGravity(double *arrCentreGrav)
{

}

// вычисление координат центра тяжести  в собственной системе координат сложного тела
void TSimpleBody_3D::calcCentreOfGravityComplicatedAxes(double *arrCentreGrav)
{
  double arrCentreGrav0 [3] = {0.};
  calcCentreOfGravity(arrCentreGrav0) ;
  mPlane.transform_xyzSKP_to_xyzSSK(arrCentreGrav0, arrCentreGrav) ;
}
// вычисление тензора моментов инерции  ля единичной массы
void TSimpleBody_3D::calcInertiaMtrx(double *arrInertTens)
{

}


// вычисление матрицы моментов инерции постого тела для массы mM
// задано простое тело  TSimpleBody_3D
// требуется вычислить матрицу моментов  инерции тела относительно
// прямоугольной декартовой системы координат, с осями параллельными осям
// начальной связанной  сиситемы координат  с центром в точке центра масс тела
// Бухгольц Основной курс теорет миеханики, ч 2, стр 136
//  делается поворот и сдвиг по теореме Штейнера
// OUTPUT:
// arrInertTens[9] - матрица моментов инерции
void TSimpleBody_3D::calcInertiaMtrxMass( double *arrInertMtrx)
{
 double arr0[9] = {0.},arr1[9] = {0.},arr2[9] = {0.},arrInertiaTensor[9] = {0.};
 // матрица моментов инерции относительно главных осей инерции  для простого телв единичной массы
 calcInertiaMtrx(arr2);
 ///

 // умножаем на массу
  MatrxMultScalar(arr2, 3, 3, mM,arr0);
  ///

  // создаем тензор инерции
  memcpy(arrInertiaTensor, arr0, 9 * sizeof(double));
  for (int i =0; i < 3; i++)
  {
	 for (int j =0; j < 3; j++)
	 {
	   if (j!=i)
	   {
		 arrInertiaTensor[3 * i + j] = -arrInertiaTensor[3 * i + j];
	   }
	 }
  }

 // поворот
 calcF_D_FTransp(mPlane.marrF, arrInertiaTensor,3,arr1 ) ;
 ///

 // сдвиг по теореме Штейнера

 calcInertiaMtrxShteinerSdvig(mM, mPlane.marrS0

 , arr1, arrInertMtrx)  ;

}
// вычисление матрицы моментов инерции отностиельно сиситемы коорлинат сдвинутой на вектор  arrS
// Бухгольц, ч 2, стр 136
// INPUT:
// VAlM -  масса тела
// arrS[3] -  вектор положения начала координат системы в которой	 надо вычислить матрицу моментов  мнерции
// arrInertTensInp[9] -  матрица ммоментов инерции в старой системе координат
// OUTPUT:
// arrInertTensOut - -  матрица ммоментов инерции в новой  системе координат
void TSimpleBody_3D::calcInertiaMtrxShteinerSdvig(const double VAlM, double *arrS
 , double *arrInertMtrxInp, double *arrInertMtrxOut)
 {

  arrInertMtrxOut[0] =  arrInertMtrxInp[0] +  VAlM*(arrS[1] * arrS[1] + arrS[2] * arrS[2]);
  arrInertMtrxOut[1] =  arrInertMtrxInp[1] +  VAlM*(arrS[0] * arrS[1]) ;
  arrInertMtrxOut[2] =  arrInertMtrxInp[2] +  VAlM*(arrS[0] * arrS[2]) ;

  arrInertMtrxOut[3] =  arrInertMtrxOut[1];
  arrInertMtrxOut[4] =  arrInertMtrxInp[4] +  VAlM*(arrS[0] * arrS[0] + arrS[2] * arrS[2]);
  arrInertMtrxOut[5] =  arrInertMtrxInp[5] +  VAlM*(arrS[1] * arrS[2]) ;

  arrInertMtrxOut[6] =  arrInertMtrxOut[2];
  arrInertMtrxOut[7] =  arrInertMtrxInp[5] ;
  arrInertMtrxOut[8] =  arrInertMtrxInp[8] +  VAlM*(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
}
#pragma package(smart_init)
