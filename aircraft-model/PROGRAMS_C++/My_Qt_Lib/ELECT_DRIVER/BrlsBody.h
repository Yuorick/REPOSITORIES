#ifndef BRLSBODY_H
#define BRLSBODY_H
#include "Brls.h"
#include "URPolygon.h"
#include "BeamPlot.h"

#define QUANT_PLOTS 4

class QBrls;
class TURPolygon;
class QBeamPlot ;
class QBrlsBody : public QBrls
{
public:
    QBeamPlot marrPlots[QUANT_PLOTS];
    int mQuantPlots;
    double marrBeginPlot[QUANT_PLOTS];
    // cellsize
    double mcs;

    // ~QBrlsBody();

    QBrlsBody();

    QBrlsBody (const QBrlsBody &R);

    QBrlsBody &operator=(const QBrlsBody  &R);


   // QBrlsBody( const double  Cv, const double   Cx, const double  VAlM
        //       , const double  VAlLengthTotal,const double  VAlHeight
         //      ,const double  Max_dOm_po_dt, const QSegment SgmBarrier
         //     , const TPlane Plane, const QBeamPlot *arrPlots
           //    ,const double* arrBeginPlot);

    QBrlsBody( const double  Cv, const double   Cx, const double  VAlM
               , const double  VAlLengthTotal,const double  VAlHeight
               ,const double  Max_dOm_po_dt, const QSegment SgmBarrier
               , const TPlane Plane, const double  cs,const int ITypeOfBRLS);

    double calcMass();

    double calcJY0();

    void   createBodyGraph(wchar_t *wchOutPutFold);

    void analysisWind(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
                      ,const double VAlStep,  double *pvalFSum
                    , double *pvalMSum,double *arrBuff, int *pquantDoneSteps);

    double calcIqu(const double VAlOmega,  const double VAlWindV
                   ,const double VAlStep, const double x);

    double calcKqu(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
                   ,const double VAlStep, const double x);

    double calcRqu(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
                   ,const double VAlStep, const double x);

    double calc_Wind_qu(const double VAlOmega, const double VAlWindV, const double VAlz);

    void createSideView(const double scale_x,const double  scaley,wchar_t *wchOutPutFold);

    void analysisIce(const double VAlPIce,const double VAlStep, double *pvalFSum
           , double *pvalMSum,double *arrBuff, int *pquantDoneSteps);

    static double calc_Hail_Mass(const double VAlDiam);

    static double calc_Hail_V(const double VAlDiam);

    static double calc_Hail_EffectH(const double VAlDiam);

    void analysisHail_(const double VAlDiam,const double VAlHailZ0,const double VAlStep, double *pvalFSum
           , double *pvalMSum,double *arrBuff, int *pquantDoneSteps);


    int calcIntervalNum(const double VAlHailZ0);

    void calcProgib(const double VAlHailZ0,const double VAlZ, const double VAlP
                               ,double *pvalProgib, double *pvalDerivProgib);

    void calcProgibInside(const double VAlHailZ0,const double VAlZCur, const double VAlP
                                ,double *pvalProgib, double *pvalDerivProgib);

};



#endif // BRLSBODY_H
