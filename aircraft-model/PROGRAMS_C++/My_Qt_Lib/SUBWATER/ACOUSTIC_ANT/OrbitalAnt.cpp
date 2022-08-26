#include "OrbitalAnt.h"
#include <math.h>
#include "Comp.h"
#include "Receiver.h"
#include "Gauss.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

using namespace std;

QOrbitalAnt::QOrbitalAnt():QSubWaterAnt()

{
  mr = 0;
  mscatter_a   = 0;
  mscatter_ang  = 0;
  mscatter_radius  = 0;
}
// конструктор копирования
 QOrbitalAnt ::QOrbitalAnt (const QOrbitalAnt &R):QSubWaterAnt(R)
 {     
     mr = R.mr;
     mscatter_a = R.mscatter_a;
     mscatter_ang = R.mscatter_ang;
     mscatter_radius  = R.mscatter_radius ;
 }


// оператор присваивания
QOrbitalAnt &QOrbitalAnt::operator=(const QOrbitalAnt  &R)
{
    if(this == &R)
    {
        return *this;
    }
    QSubWaterAnt:: operator= (R);
    mr = R.mr;
    mscatter_a = R.mscatter_a;
    mscatter_ang = R.mscatter_ang;
    mscatter_radius  = R.mscatter_radius ;
  return *this ;
}
// парам констр
QOrbitalAnt :: QOrbitalAnt( const QVector<QReceiver> vctOrbitalAnt, const QVector<QPointF> vctPntXY)

{
   mvctReceiver = vctOrbitalAnt ;
   mvctPntXY = vctPntXY;
}

//----------------------------------------
// парам констр
 QOrbitalAnt :: QOrbitalAnt(const int QuantReceivers, const double a_Receiver
                            , const double radius):QSubWaterAnt (QuantReceivers,  a_Receiver
                                                                 ,  radius)

{
     mr = radius;
fillReceiversLocation( radius);

}

//--------------------------------
void QOrbitalAnt ::fillReceiversLocation(const double radius)
{
    double step = 2. * M_PI/ ((double)mQuant);
    for (int i = 0; i < mQuant; ++i)
    {
        double ang = ((double)i) * step;
        mvctPntXY.replace(i, QPointF(radius * sin(ang),radius * cos(ang)));
    }
}


