#ifndef OrbitalAnt_H
#define OrbitalAnt_H
#include <QVector>
#include <QPointF>
#include "SubWaterAnt.h"
class QReceiver;
class TComp;


class QOrbitalAnt:public QSubWaterAnt
{
public:
    QOrbitalAnt();

    QOrbitalAnt  (const QOrbitalAnt &R) ;
    // оператор присваивания
    QOrbitalAnt &operator=(const QOrbitalAnt  &R);
    // парам констр
    QOrbitalAnt( const QVector<QReceiver> vctOrbitalAnt, const QVector<QPointF> vctPntXY);

    QOrbitalAnt(const int QuantReceivers, const double a_Receiver
                                      , const double radius);





    double mr;
    double mscatter_a;
    double mscatter_ang;
    double mscatter_radius;

    virtual void fillReceiversLocation(const double radius);


};

#endif // OrbitalAnt_H
