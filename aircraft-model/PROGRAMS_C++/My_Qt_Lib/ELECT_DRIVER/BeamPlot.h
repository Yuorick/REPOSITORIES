#ifndef BEAMPLOT_H
#define BEAMPLOT_H
#include "URPolygon.h"

#define MAX_QUANT_PROFILES 4

class TURPolygon;
class QBeamPlot
{
public:
    TURPolygon marrPlgProfile[MAX_QUANT_PROFILES];
    double marrRo[MAX_QUANT_PROFILES];
    double marrSigma[MAX_QUANT_PROFILES];
    double marrE[MAX_QUANT_PROFILES];
    int mquantPlg;
    double marrPseudoCentre[2];
    double marrMtrxPseudoInertia[4];
    double mcs;

  //  ~QBeamPlot();

    QBeamPlot();

    QBeamPlot (const  QBeamPlot &R);

    QBeamPlot  &operator=( const QBeamPlot  &R);

    QBeamPlot (const TURPolygon Profile , const double  Ro, const double  s);

    QBeamPlot &operator=(const TURPolygon  &R);

    QBeamPlot(  const int quantPlg, const TURPolygon *arrPlgProfile, const double *arrRo
                , const double *arrSigma, const double *arrE, const double mcs);

  //  QBeamPlot(const int INumPlotFurke2540);

   static void createBeamPlotFurke2540(const double cs , const int INumPlotFurke2540, QBeamPlot &BeamPlot);
//----------------------------------------------------------------------

    double calcNeutral_Y();

    double calcNeutral_X();

    void   calcPseudoInertiaMtrx( double *arrPseudoCentre, double *mtrxPseudoInertia);

    double calcCurveRad_PlaneZX(const double VAlMomY);

    double calcCurveRad_PlaneZY(const double VAlMomX);

    double cal_LinDensityMass();

    double calJY0(const double a, const double b);

    void createGraph(wchar_t *wchFoldCur);

    void calcCoeffMaxSigma_ZX( double *arrSigmaCoeff );

    static double max__(const double a, const double b);

    void createSideView(const double a,const double  b
                        ,const double scale_x,const double  scaley,wchar_t *wchFileCur);

    void calcCoeffMaxSigma_ZY( double *arrSigmaCoeff );





};

#endif // BEAMPLOT_H
