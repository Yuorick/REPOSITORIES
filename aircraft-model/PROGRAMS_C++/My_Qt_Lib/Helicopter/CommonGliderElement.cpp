#include "CommonGliderElement.h"
#include <string.h>


TCommonGliderElement::TCommonGliderElement()
{

}

TCommonGliderElement ::TCommonGliderElement(const long double S
                                          , const TLongPlane PlaneSvSK, const long  double Cx0, const long  double Cy0, const long double AlfCrit
                                                                                   ,long  double *arrAirP)
    :TAbstractGliderElement( S, PlaneSvSK, Cx0, Cy0,  AlfCrit , arrAirP)
{

}

TCommonGliderElement ::TCommonGliderElement(long double *arrInpDataGliders): TAbstractGliderElement(arrInpDataGliders)
{

}

// конструктор копирования
 TCommonGliderElement ::TCommonGliderElement (const TCommonGliderElement &R)
 {
     mS = R.mS;
     mPlaneSvSK   = R.mPlaneSvSK ;
     mCx0  = R.mCx0 ;
     mCy0 = R.mCy0;
     mAlfCrit = R.mAlfCrit;
     memcpy(marrAirP, R.marrAirP,  3 * sizeof(long double));
 }
 // оператор присваивания
 TCommonGliderElement TCommonGliderElement::operator=(TCommonGliderElement  R)
 {
     mS = R.mS;
     mPlaneSvSK   = R.mPlaneSvSK ;
     mCx0  = R.mCx0 ;
     mCy0 = R.mCy0;
     mAlfCrit = R.mAlfCrit;
     memcpy(marrAirP, R.marrAirP,  3 * sizeof(long double));
    return *this ;
  }
// вычисление матрицы разворота СвСК элемента планера
// для элемента типа руля эта матрица считакется,
// для простого элемента приравнивается единичной
// INPUT:
// VAlAlfa - угол поворота руля
//OUTPUT:
// arrMwave[9] - матрица 3х3
void TCommonGliderElement::calMtrxMwave(const long  double VAlAlfa,long  double *arrMwave)
{
 memset(arrMwave, 0, 9 * sizeof(long double));
 arrMwave[0] = arrMwave [4]= arrMwave[8] = 1.;

}
