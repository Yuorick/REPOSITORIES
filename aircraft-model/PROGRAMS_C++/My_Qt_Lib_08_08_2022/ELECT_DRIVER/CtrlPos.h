#ifndef CTRLPOS_H
#define CTRLPOS_H
#include "Ctrl.h"
#include "CtrlFollow.h"
#include "MeasurmentImitator.h"

class QMeasurmentImitator;

class QCtrlPos:public QCtrlFollow
{
public:
    QCtrlPos();
    QCtrlPos(const QCtrlPos &R);
    QCtrlPos  &operator=( const QCtrlPos  &R);

    QCtrlPos( const QElectMotor ElectMotor ,  QLoad*  Load, const double *arrSpreadParams
              , const double MomOut, const double VAlTettaBegin,const double VAlOmegaBegin
              , const double T0,const double TCur
              ,const double valh,TComp *pCmpArrLamb);



};

#endif // CTRLPOS_H
