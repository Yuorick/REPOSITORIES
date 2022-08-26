#ifndef DRIVEMOVEIMIT_H
#define DRIVEMOVEIMIT_H
#include "Ctrl.h"
#include "ElectDriver.h"
#include "Environment.h"
#include "MeasurmentImitator.h"
#include "Comp.h"
#include "AbstractTraj.h"

class QCtrl;
class QElectDriver;
class TEnvironment;
class mMeasurmentImitator;
class TComp;
class QAbstractTraj;

class QDriveMoveImit
{
public:
    QCtrl *mpCtrl;
    QElectDriver mRealDriver;
    QMeasurmentImitator mMeasurmentImitator;
    // ветер
    TEnvironment mEnvironment;

    // траектория цели
    QAbstractTraj *mpTargTraj;

    // внешний момент реальный
    double mMomOutReal;

    // внешний момент модельный
    double mMomOutModel;



    QDriveMoveImit();

    QDriveMoveImit (const  QDriveMoveImit &R);

    QDriveMoveImit  &operator=( const QDriveMoveImit &R);

    QDriveMoveImit ( QCtrl *pCtrl, const QElectDriver RealDriver,
      const QMeasurmentImitator MeasurmentImitator,  QAbstractTraj *pTargTraj
                     ,const double MomOutReal,const double MomOutModel);

    QDriveMoveImit ( QCtrl *pCtrl, const QElectDriver RealDriver,
          const QMeasurmentImitator MeasurmentImitator, const TEnvironment Environment
        ,  QAbstractTraj  *pTargTraj,const double MomOutReal,const double MomOutModel);
    //-----------------------------------------------------------
  //  void move(const double VAlTBegin,const double VAlMovingT,const double VAlIntegrStep0
     //         ,double *arrBuff, int *piQuantRows);

    void move(double *arrQObjective,const double VAlTBegin,const double VAlMovingT,const double VAlIntegrStep0
               ,double *arrBuff, int *piQuantRows);

    void StatisticProcess(double *arrBuff, const int quantDoneSteps,const double VAl_TPeriod, double *arrE, double *arrDisp);

    void calcAnaliticSumTrajScatters(double *arrObjective,double *arrW,double *arr_dx_po_dAlf, double *arrAHarm
                                     , double *arrSystErr,double *arrERandomSig, double *arrSumSig);
    void calcSystTrajScatters(double *arrObjective,double *arrW,double *arr_dx_po_dAlf,  double *arrSystErr) ;

    void createVectParamsDiffers(double *arrParamsDiffer) ;

    void calcHarmAmpVect(double *arrObjective,double *arrW, double *arrAHarmSquare );

    void calcFluctK(double *arrW,  double *arrFluktK);


};

#endif // DRIVEMOVEIMIT_H
