#ifndef LINEXY_H
#define LINEXY_H
#include <QLineF>
class QLineF;
class QPointF;
class QPointXY;
class QRectF;
class TURPolyLine;

class QLineXY : public QLineF
{
public:

    QLineXY();

    QLineXY (const  QLineXY &R);

    QLineXY  &operator=( const QLineXY  &R);

    QLineXY(const QPointF &p1, const QPointF &p2);

    QLineXY(qreal x1, qreal y1, qreal x2, qreal y2);

    QLineXY(const QLine &line);

    QPointXY getPXY1();

     QPointXY getPXY2();

    static void createVectArrow(const double xbegin, const double ybegin
                                , const double e_directx,  const double e_directy
                                ,const double length,const double arrowLength
                                , const double valAng, QVector<QLineF> *vctArrow);

    static double Sign(const double x);

    static void rastiagenie(QVector<QLineF> *vctLine, const double coef);

    void rastiag( const double coef);

    static void translateVct(QVector <QLineF> *vctLine, QPointF &pntOffset);

    void stretch(const double &scalx,const double &scaly);

    void stretch(QVector<QLineF> *vctLine,const double &scalx,const double &scaly);

    static void stretch(QLineF& ln0,const double &scalx,const double &scaly);

    static void stretchVct(QVector<QLineF> *vctLine,const double &scalx,const double &scaly);

    static void createVectNet(const double xmin, const double xmax
                                  ,const double ymin, const double ymax
                                  , const double xstep,  const double ystep
                                  , QVector<QLineF> *vctNet);

    static void createVectAxes(const double xmin, const double xmax
                                  ,const double ymin, const double ymax
                                 ,const double arrowLength, const double valAng
                                  , QVector<QLineF> *VctAxes);

    static QRectF  boundBox(QLineF &ln);

    bool   isIntersected(const QRectF &Rect);

    static bool   isIntersected(const QVector<QLineF> &vctLine, QRectF &Rect);

    static QLineXY createLineXY(const QLineF &ln);

    static TURPolyLine createTURPln(const QLineF &ln);

    static TURPolyLine createTURPln(const QVector<QLineF> &vctLine);

};

#endif // LINEXY_H
