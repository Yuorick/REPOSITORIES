#ifndef CTRLVELOQUIQ_H
#define CTRLVELOQUIQ_H
#include "CtrlVelo.h"
#include "MeasurmentImitator.h"

class QCtrlVelo;
class QVeloQuiqData;
class QMeasurmentImitator ;

class QVeloQuiqData
{
    public:
    // момент переключения
    double mOverSwitchT;

    // конечное время
    double mEndT;

    //
    double marrObjective[2];

    int misignum;
     QVeloQuiqData();

     QVeloQuiqData(const QVeloQuiqData &R);

     QVeloQuiqData  &operator=( const QVeloQuiqData  &R);


};

class QCtrlVeloQuiq: public QCtrlVelo
{
public:

    // ограничение на максимум напряжения Uq
    double mMaxU;
    QVeloQuiqData mVeloQuiqData;


    QCtrlVeloQuiq();

    QCtrlVeloQuiq(const QCtrlVeloQuiq &R);

    QCtrlVeloQuiq  &operator=( const QCtrlVeloQuiq  &R);

    QCtrlVeloQuiq( const QElectMotor ElectMotor ,  QLoad*  Load
                 , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                 ,const double VAlOmegaStatBegin , const double T0,const double TCur
                 , const double *arrObjective
                 ,const double MaxU,const double valh,TComp *pCmpArrLamb);

    virtual void calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU);

    bool calcSwithOverTimes(QVeloQuiqData *pVeloQuiqData,const double VAlMomOut);
};


#endif // CTRLVELOQUIQ_H
