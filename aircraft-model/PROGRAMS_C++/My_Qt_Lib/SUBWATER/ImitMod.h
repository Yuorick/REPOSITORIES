#ifndef IMITMOD_H
#define IMITMOD_H
#include "PeaceVess.h"
#include "Table_1D.h"

class QPeaceVess;
class QBigMeasure;
enum enumTypeOfVessTraj{ PNT_6, PNT_7,QDRT, ZIG_ZAG, DIAMETRS, LINE6};


class QImitMod
{
public:

    QImitMod();

    QImitMod (const  QImitMod &R);

    QImitMod  &operator=( const QImitMod  &R);

    // парам конструктор 1
    QImitMod (const QPeaceVess Vess, const TTable_1D tblRealPrfl
              ,const double *arrTrueBeaconPos,const double *arrAprioriBeaconPos
              , const enumTypeOfVessTraj TypeOfVessTraj, const double TObrabotki
                        ,const double *arrAprioriGpsPosParams);

    // парам конструктор 2
    QImitMod (const QPeaceVess Vess, const TTable_1D tblRealPrfl
              ,const double *arrTrueBeaconPos,const double *arrAprioriBeaconPos
              , const enumTypeOfVessTraj TypeOfVessTraj, const double TObrabotki
           ,const double *arrAprioriGpsPosParams, const double *arrTrueBeaconVelo, const double *arrAprioriBeaconVelo );

    void createMeasures (const int QuantMeas, QBigMeasure *p);

    bool createSingleMeasure (  QBigMeasure &Meas);

    void imitateMeasure( QTrueMeasParams &trueMeasParams, QBigMeasure &Meas);

    void imitateVesselPosition(double *arrSGskVessTrue,double *arrMu,double *arrMuZv,double *arrSGskVessZv);

    void  initPartOfQdrtTraj(const int numQdrtSide
                             ,const double *arrAprioriBeaconPos,const double ZonaR);

    int calcNumQdrtSide(const double *arrQdrtCentre, const double VAlQdrtSideLength
                                  , const int NUmQdrtSideCur);

    void createSeparatedPointsMeasures( const int QuantMeas, const double ZonaR_
                                        , QBigMeasure *arrMeas);

    void createQdrtTrajMeasures( const int QuantMeas, const double ZonaR_
                            , QBigMeasure *arrMeas);

    void createZigZagTrajMeasures( const int QuantMeas, const double ZonaR_
                            , QBigMeasure *arrMeas);

    int calcNumPartOfZigZagTraj(const double *arrQdrtCentre, const double VAlQdrtSideLength
                                  , const int NUmQdrtSideCur);

    void  initPartOfZigZagTraj(const int numQdrtSide
                    ,const double *arrAprioriBeaconPos,const double VAlQdrtSideLength);

    void createDiamTrajMeasures( const int QuantMeas, const double ZonaR_
                            , QBigMeasure *arrMeas);

    void  initPartOfDiamTraj(const int numDiamPart
                    ,const double *arrAprioriBeaconPos,const double VAlDiamLength);

    void createLine6TrajMeasures( const int QuantMeas, const double ZonaR_
                           , QBigMeasure *arrMeas);

    void  initPartOfLine6Traj(const int numDiamPart
                    ,const double *arrAprioriBeaconPos,const double VAlDiamLength);

    QPeaceVess mVess;

    TTable_1D mtblRealPrfl;

    // координаты маяка  истинные
   double marrTrueBeaconPos[3];
   // координаты маяка  априорные
   double marrAprioriBeaconPos[3];

   // скрость маяка истиная
  double marrTrueBeaconVelo[3];
  // координаты маяка  априорные
  double marrAprioriBeaconVelo[3];

    enumTypeOfVessTraj mTypeOfVessTraj;

    double mTObrabotki;

    // вектор параллакса GPS априорный
    double marrAprioriGpsPosParams[6];



};

#endif // IMITMOD_H
