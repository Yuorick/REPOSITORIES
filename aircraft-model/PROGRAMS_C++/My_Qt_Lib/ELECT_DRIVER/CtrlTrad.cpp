#include "CtrlTrad.h"
#include <math.h>
#include <string.h>

QCtrlTrad::QCtrlTrad():QCtrl()
{
 mModU = 0.;
 mT1 = 0.;
 mT2 = 0.;
}

// конструктор копирования
 QCtrlTrad::QCtrlTrad(const QCtrlTrad &R):QCtrl( R)
 {
  mModU = R.mModU;
  mT1 = R.mT1;
  mT2 = R.mT2;
 }
// оператор присваивания
 QCtrlTrad  &QCtrlTrad::operator=( const QCtrlTrad  &R)
{
     if(this == &R)
     {
         return *this;
     }
    mModU = R.mModU;
    mT1 = R.mT1;
    mT2 = R.mT2;

    QCtrl:: operator= (R);

    return *this ;
}
//-----------------------------------

// парам конструктор
 // параметрический конструктор
QCtrlTrad::QCtrlTrad( const QElectMotor ElectMotor ,  QLoad*  Load
                          , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                          ,const double VAlOmegaStatBegin , const double T0,const double TCur
                      ,const double valh,TComp *pCmpArrLamb, const double ModU, const double T1, const double T2)
     :QCtrl (ElectMotor ,  Load, arrSpreadParams, MomOut, VAlTettaBegin
              ,VAlOmegaStatBegin ,T0,  TCur,valh,pCmpArrLamb)
 {
   mModU = ModU;
   mT1 = T1;
   mT2 = T2;
 }

//------------------------------------------------------------------
void QCtrlTrad::calcCurU(const double *arrObjective,const double VAlTObjective, double *arrU)
{
    double arrEstPhVect[QVARS] = {0.};
    memcpy(arrEstPhVect, mFiltr.marrCurEst, sizeof(double) * QVARS );
    double arrStatorU[2] = {0.};
    double valModU = mModU;
    if ((mTCur - mT0) > mT1)
    {
        valModU = 0.;
    }
    arrStatorU[0] = valModU * cos((mTCur - mT0) * arrObjective[1] );
    arrStatorU[1] = valModU * sin((mTCur - mT0) * arrObjective[1] );

    arrU[0] = cos(arrEstPhVect[3]) * arrStatorU[0] + sin(arrEstPhVect[3]) * arrStatorU[1] ;
    arrU[1] = -sin(arrEstPhVect[3]) * arrStatorU[0] + cos(arrEstPhVect[3]) * arrStatorU[1] ;


}
