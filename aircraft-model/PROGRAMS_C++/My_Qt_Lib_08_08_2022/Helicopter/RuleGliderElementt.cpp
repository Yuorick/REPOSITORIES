#include "RuleGliderElement.h"
#include <string.h>
#include <math.h>


TRuleGliderElement::TRuleGliderElement()
{

}

TRuleGliderElement ::TRuleGliderElement(const long double S
                      , const TLongPlane PlaneSvSK, const long  double Cx0, const long double Cy0, const long double AlfCrit
                      ,long double *arrAirP)
                      :TAbstractGliderElement( S, PlaneSvSK, Cx0, Cy0,  AlfCrit , arrAirP)
{

}

TRuleGliderElement ::TRuleGliderElement(long double *arrInpDataGliders): TAbstractGliderElement(arrInpDataGliders)
{

}

// конструктор копирования
 TRuleGliderElement ::TRuleGliderElement (const TRuleGliderElement &R)
 {
     mS = R.mS;
     mPlaneSvSK   = R.mPlaneSvSK ;
     mCx0  = R.mCx0 ;
     mCy0 = R.mCy0;
     mAlfCrit = R.mAlfCrit;
     memcpy(marrAirP, R.marrAirP,  3 * sizeof(long double));
 }
 // оператор присваивания
 TRuleGliderElement TRuleGliderElement::operator=(TRuleGliderElement  R)
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
 // это матрица смостоящая из столбцов един. ортов развернутой СвСК
 // в осях начальной СвСК
// INPUT:
// VAlAlfa - угол поворота руля
//OUTPUT:
// arrMwave[9] - матрица 3х3
void TRuleGliderElement::calMtrxMwave(const long double VAl_delta, long double *arrMwave)
{
 arrMwave[0] = cosl(VAl_delta);
 arrMwave[1] = -sinl(VAl_delta);
 arrMwave[2] = 0.;

 arrMwave[3] = sinl(VAl_delta);
 arrMwave[4] = cosl(VAl_delta);
 arrMwave[5] = 0.;

 arrMwave[6] = 0.;
 arrMwave[7] = 0.;
 arrMwave[8] = 1.;
}
