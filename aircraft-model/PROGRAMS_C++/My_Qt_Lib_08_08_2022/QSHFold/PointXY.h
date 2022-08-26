#ifndef POINTXY_H
#define POINTXY_H
#include <QPointF>


class QPointXY : public QPointF
{
public:
    QPointXY();

    QPointXY(const QPoint &point);

    QPointXY(qreal xpos, qreal ypos);

    QPointXY (const QPointXY &R2);

    QPointXY &operator=(const QPointXY  &R2);

    QPointXY   LinTransform(const double  valAng , const QPointXY pntSdvig,const double valRastigenie );

    static QPointXY toPntXY(const QPointF& pntF);

    static double calcVectS(const QPointF& P1,const QPointF& P2);

    static double calc_dist(const QPointF& P1, const QPointF& P2);
};

#endif // POINTXY_H
