#ifndef PARAMSID_H
#define PARAMSID_H
#include "ElectDriver.h"
#include "MeasurmentImitator.h"


class QMeasurmentImitator;
class QElectDriver;
class QStatSolutionParams;
class QRezPointTraj;

class QParamsID
{
public:
    QElectDriver mDriverModel;

    QElectDriver  mDriverReal;

    QMeasurmentImitator mMeasurmentImitator;

    double marrSpreadParams[QPARAMS];


    // массив номеров идентифицируемых параметров
    int miarrTargNums[QPARAMS + 2];
    // к-во идентифицируемых параметров
    int mlenTargNums;

    QParamsID();

    QParamsID (const QParamsID &R);

    QParamsID &operator=(const QParamsID  &R);

    QParamsID(const QElectDriver DriverModel,const QElectDriver DriverReal
                , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
                , const int *iarrTargNums , const double lenTargNums );


    static QElectDriver static_identify_by_StatMeth(const QElectDriver DriverModel,const QElectDriver DriverReal
                                             , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
                                             ,const double DEltaTetta, const int QUantTetta,const double VAlOmega0,const double DEltaOmega, const int QUantOm
                                             ,const double  VAlMovingT,TComp *cmpArrLamb, const int QUantLambSets,int *iarrTargNums
                                             ,const int LEnTargNums, const double VAlIntegrStep );

    void imitateStatExperimentsRez(const double DEltaTetta, const int QUantTetta,const double DEltaOmega,const double VAlOmega0
                                   , const int QUantOm,const double  VAlMovingT, TComp *cmpArrLamb, const int QUantLambSets
                                     , QStatSolutionParams *parrStatSolutionParams,double *arrDynamicDelX, const double VAlIntegrStep );

    double calcDelAlf_StatMeth(int quantExprs, QStatSolutionParams *parrStatSolutionParams
                               , double *arrDynamicDelX,  double *arrDElAlf);

    static QElectDriver identify_by_DynamicMeth(const QElectDriver DriverModel,const QElectDriver DriverReal
                                                                       , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
                                                                       ,const double  VAlU0,const double   VAlT1,const double VAlMovingT,const double TrgOm, const double TetRotor0
                                                                       , const int QuantIsp, int *iarrTargNums,const int LEnTargNums, const double VAlIntegrStep,double *arrBuff
                                                                       , int *pquantDoneSteps);

    void imitateDynamicExperimentsRez(const double  VAlU0, const double VAlT1
                                      ,const double  VAlMovingT,const double TrgOm, const double TetRotor0
                                      , const int QUantIsp,const double  VAlIntegrStep ,double *arrBuff,int *quantDoneSteps);

    QElectDriver identify_by_DynamicMeth_( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
            ,const int quantDoneSteps);

    void identify_by_DynamicMeth_1( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
                                    ,const int quantDoneSteps, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut);

    double identify_by_DynamicMeth_2( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
                                    ,const int quantDoneSteps, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut);

    double identify_by_StatMeth(int quantExprs, QStatSolutionParams *parrStatSolutionParams
                                , double *arrDynamicDelX, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut);
};


class QRezPointTraj
{
public:
  double mt; // время привязки точки траектории
  double marrX[QVARS];// фазовый вектор
  double marrU[2];  // вектор управлений

  QRezPointTraj();
  QRezPointTraj (const QRezPointTraj &R);
  QRezPointTraj &operator=(const QRezPointTraj  &R);
  QRezPointTraj(  const double t, const double *arrX ,  const double *arrU );
};

#endif // PARAMSID_H
