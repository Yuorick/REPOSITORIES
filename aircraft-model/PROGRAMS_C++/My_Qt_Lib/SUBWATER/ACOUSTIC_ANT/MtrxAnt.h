#ifndef MtrxAnt_H
#define MtrxAnt_H
#include <QVector>
#include <QPointF>
#include "SubWaterAnt.h"
class QReceiver;
class TComp;


class QMtrxAnt:public QSubWaterAnt
{
public:
    QMtrxAnt();

    QMtrxAnt  (const QMtrxAnt &R) ;
    // оператор присваивания
    QMtrxAnt &operator=(const QMtrxAnt  &R);
    // парам констр
    QMtrxAnt( const QVector<QReceiver> vctMtrxAnt, const QVector<QPointF> vctPntXY);

    QMtrxAnt(const int QuantReceivers, const double a_Receiver
                                      , const double radius);





    double mr;
    double mscatter_a;
    double mscatter_ang;
    double mscatter_radius;

    virtual void fillReceiversLocation(const double radius);


};

#endif // MtrxAnt_H
