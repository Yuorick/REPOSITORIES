#ifndef CTRLTRAD_H
#define CTRLTRAD_H
#include "Ctrl.h"
#include "MeasurmentImitator.h"

class QCtrlTrad:public QCtrl
{
public:
    QCtrlTrad();

    QCtrlTrad(const QCtrlTrad &R);

    QCtrlTrad  &operator=( const QCtrlTrad  &R);

    QCtrlTrad( const QElectMotor ElectMotor ,  QLoad*  Load
               , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
               ,const double VAlOmegaStatBegin , const double T0,const double TCur
           ,const double valh,TComp *pCmpArrLamb, const double ModU, const double T1, const double T2);

    // модуль напряжкения
    double mModU;
    // интервал времени действия напряжения [0; T1]
    double mT1;
    //интервал времени нулевого напряжения  [T1;T2]
    double mT2;

    virtual void calcCurU(const double *arrObjective,const double VAlTObjective, double *arrU);
};

#endif // CTRLTRAD_H
