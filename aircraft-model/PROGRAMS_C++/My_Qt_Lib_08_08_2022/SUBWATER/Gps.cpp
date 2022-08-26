#include "Gps.h"
#include  <string.h>
#include "Gauss.h"

//
 QGps ::QGps ()
 {
     mSigXY = 0.;
     memset(marrPos, 0, 3 * sizeof(double));
 }
 // конструктор копирования
  QGps ::QGps (const QGps &R)
  {
      mSigXY = R.mSigXY ;
      memcpy(marrPos, R.marrPos, 3* sizeof(double));
  }


 // оператор присваивания
 QGps &QGps::operator=(const QGps  &R)
 {
   mSigXY = R.mSigXY ;
   memcpy(marrPos, R.marrPos, 3 * sizeof(double));
   return *this ;
 }
// парам констр
QGps :: QGps(  const double SigXY,const double* arrPos)

 {
    mSigXY = SigXY ;
    memcpy(marrPos, arrPos, 3 * sizeof(double));
}
//--------------------------------------
// имитация замера
void QGps::imitateMeasure(double *arrTruePosition_GSK, double *arrPosition_GSK_Zv)
{
     // ТЕСТ
    for (int i =0; i < 3; ++i)
    {
      arrPosition_GSK_Zv[i] =  arrTruePosition_GSK[i] +  getGauss(0, mSigXY );
    }
}
