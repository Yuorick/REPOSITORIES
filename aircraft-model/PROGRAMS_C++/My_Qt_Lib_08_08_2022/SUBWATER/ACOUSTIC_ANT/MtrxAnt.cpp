#include "MtrxAnt.h"
#include <math.h>
#include "Comp.h"
#include "Receiver.h"
#include "Gauss.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

using namespace std;

QMtrxAnt::QMtrxAnt():QSubWaterAnt()

{
  mr = 0;
  mscatter_a   = 0;
  mscatter_ang  = 0;
  mscatter_radius  = 0;
}
// конструктор копирования
 QMtrxAnt ::QMtrxAnt (const QMtrxAnt &R):QSubWaterAnt(R)
 {     
     mr = R.mr;
     mscatter_a = R.mscatter_a;
     mscatter_ang = R.mscatter_ang;
     mscatter_radius  = R.mscatter_radius ;
 }


// оператор присваивания
QMtrxAnt &QMtrxAnt::operator=(const QMtrxAnt  &R)
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
QMtrxAnt :: QMtrxAnt( const QVector<QReceiver> vctMtrxAnt, const QVector<QPointF> vctPntXY)
    :QSubWaterAnt (vctMtrxAnt,  vctPntXY)

{

}

//----------------------------------------
// парам констр
 QMtrxAnt :: QMtrxAnt(const int QuantReceivers, const double a_Receiver
                                  , const double radius):QSubWaterAnt (QuantReceivers,  a_Receiver
                                                                       ,  radius)

{
     mr = radius;
   fillReceiversLocation( radius);

}


//--------------------------------
void QMtrxAnt ::fillReceiversLocation(const double step)
{
    const int n =(int)sqrt((double)mQuant+0.01);
    double x0 = 0., y0 = 0.;
    int m = n/2;
    switch(n%2)
    {
    case 0:
        x0 = -((double)m )* step;
        y0 = ((double)m )* step;
        break;

    case 1:
        x0 = -((double)m )* step - step/2.;
        y0 = ((double)m )* step + step/2.;
        break;



     }

    for (int i = 0; i < n; ++i)
        for (int j = 0 ;j < n; ++j)
    {

        mvctPntXY.replace(i* n + j, QPointF(x0+ i * step,y0 - j * step));
    }
}


