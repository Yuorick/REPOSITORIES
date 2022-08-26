#ifndef SEGMENT_H
#define SEGMENT_H
#include "URPolyLine.h"
#include "URFigure.h"

class TURPolyLine;

class QSegment:public TURPolyLine
{
public:
    QSegment();

    QSegment(const  TURPolyLine pln);

    QSegment(const TURPointXY pnt0,const TURPointXY pnt1);

    static void calcOpen_and_Shadow_Segments( QSegment segmBase,  QSegment segmBarrier, double *arrV
                                              , QSegment *arrOpenSegm, int *pQuantOpenSegm, QSegment *pShadowSegm, int *pQuantShadowSegm );

    bool   IsVertical();

    bool  calcTang(double *pvalTang);

    bool   calcY(const double x, double *py);

    bool   calcPointXY(const double x, TURPointXY *pPointXY);

    double calcSlopeAng();


};

#endif // SEGMENT_H
