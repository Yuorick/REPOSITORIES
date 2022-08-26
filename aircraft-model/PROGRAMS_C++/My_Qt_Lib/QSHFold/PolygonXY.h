#ifndef POLYGONXY_H
#define POLYGONXY_H
#include <QLineF>
#include <QPointF>
#include "LineXY.h"
#include "PointXY.h"
#include <QPolygonF>

class TURPolygon;
class TURPolyLine;

class QPolygonXY : public QPolygonF
{
public:
    QPolygonXY();
    QPolygonXY(int size);
    QPolygonXY(const QVector<QPointF> &points);
    QPolygonXY(QVector<QPointF> &&v);
    QPolygonXY(const QRectF &rectangle);
    QPolygonXY(const QPolygon &polygon);
    QPolygonXY(const QPolygonF &polygon);
    QPolygonXY(QPolygonF &&other);

    double Box[4];
    //-----------------------------------
    //-----------------------------------
    void calcBoundBox();

    QPolygonXY(const QVector<QPointXY> &points);

    QPolygonXY(const QVector<QPointXY> points);



    QPolygonF createPlgF();

    static QPolygonXY createPlgXY(const QPolygonF &plgF);

    QPolygonXY   LinTransform(const double  valAng , const QPointXY pntSdvig
                                          ,const double valRastigenie );

    void   LinTransform_(const double  valAng , const QPointXY pntSdvig
                                          ,const double valRastigenie );

    TURPolygon createTURPlg();

    TURPolyLine createTURPln();

    QPolygonXY createPlnLowFront();

    QPolygonXY createPlnUpFront();

    double calcSq();

    bool isClockwise();

    void flipPoints();

    int findNumberUtmostRightPoint();

    static bool calcDistBetweenFrontLines(QPolygonXY& plnUp,QPolygonXY& plnLow, double &valDist);

    int findNumberUtmostLeftPoint();

    static bool calcSdvig(QPolygonXY& plnUp,QPolygonXY& plnLow, double &valSdvig);

    bool calcDistNearestRightPoint(const double x, double &valDistCur) ;

    bool calcVertDistFromPoint(const QPointF &pnt, double &valDist);

    bool findFirstIntersectPoint(const QPointF &pntNose
         ,const double* arrVIzd_TargSvSK0, QPointF &pntOutput);

    void stretch(const double &scalx,const double &scaly);

    static void stretch(QVector<QPolygonF> *vctLine,const double &scalx,const double &scaly);

    static void stretch(QPolygonF &plg,const double &scalx,const double &scaly);

    int findNumberUtmostUpPoint();

    int findNumberUtmostLowPoint();

    static void mirrowRect_X(QRectF &rect);

    static int findMostDistantVertexNumber(QPolygonF &plg,QPointF &pnt);

    bool findUtmostLowPoint_Bounded_Rect(const QRectF &Rect, QPointF *intersectionPoint);

    void  findIntersectingRectPointsArray(const QRectF &Rect, QVector <QPointF> *vctIntersect);

    void  findIntersectingLinePointsArray(const QLineF &Rect, QVector <QPointF> *vctIntersect);

    bool   isIntersected(const QRectF &rect) ;

    static TURPolygon createTURPlg(QVector <QPolygonF> &vctPlg);

    static void createPoints3DModel(const QPolygonF &plg,
         const  double &lenghtStep,const  double  &angleStep
        , bool &bNoseOnly,QVector <double> &buffVectPoints);

    static TURPolygon createTURPlg(const QPolygonF &plg);

    static void   flip(QPolygonF &plg);

    static double calcSq(const QPolygonF &plg);

    static bool findFirstIntersectPoint(const QPolygonF &plg,const QPointF &pntNose
                              ,const double* arrVIzd_TargSvSK0, QPointF &pntOutput);
};

#endif // POLYGONXY_H
