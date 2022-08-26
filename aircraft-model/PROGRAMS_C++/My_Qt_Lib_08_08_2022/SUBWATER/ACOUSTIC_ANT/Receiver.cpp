#include "Receiver.h"
#include <math.h>

QReceiver::QReceiver():
ma(0.014)
{

}
// конструктор копирования
 QReceiver ::QReceiver (const QReceiver &R)
 {

     ma = R.ma;
 }


// оператор присваивания
QReceiver &QReceiver::operator=(const QReceiver  &R)
{

  ma = R.ma;
  return *this ;
}
// парам констр
QReceiver :: QReceiver(  const double a)

{  
   ma = a;
}
//-----------------------------------------
// Диаграмма направленности
double QReceiver::fncBessDiagr(const double tetta, const double lamb)
{
    double valk = 2. * M_PI/ lamb;
    double x = valk * ma* sin(tetta);
    return 2. * fncBx(x);
}


//-----------------------------------------
// Диаграмма направленности
//double QReceiver::clcAmpBessDiagr(const double tetta,const double dist
    //                              ,const double transmitterWaveAmp, const double lamb)
//{
  //  double diagr = fncBessDiagr( tetta, lamb);
  //  double amp = calcSignalAmp(dist
  //                              ,transmitterWaveAmp, lamb);

  //  return fabs(amp * diagr);
//}
//----------------------------------------
// Амплитуда водны на на приемнике
double QReceiver::calcSignalAmp(const double d,const double transmitterWaveAmp, const double lamb)
{
    return transmitterWaveAmp /d * exp(-log(10)
                   * 0.036 * 30. * sqrt(30.) * d/20./1000.) ;
}




//------------------------------
// Вычисление функции Бесселя 1
double fncBessel1(const double x)
{
   double valSign= 1.;
    double sum = 1.;
    double y =  x* x /4.;
    double prod = 1.;
    for (int i = 1; i < 50; ++i)
    {
       prod = -prod *  y/ ((double) (i * (i +1))) ;
       sum += prod;
    }
    sum = sum * x /2.;

    return sum;
}
//------------------------------
// Вычисление производной  функции Бесселя1
double fnc_dBessel1_po_dx(const double x)
{
    double valSign= 1.;
    double sum = 1./2.;
    double y = x* x /4.;
    double prod = 1./2.;
    for (int i = 1; i < 50; ++i)
    {

       prod = -prod *  y/ ((double) (i * (i +1))) ;
       sum += prod * ((double)( 2 * i + 1));
    }
    return sum;
}
//------------------------
// фукнкция BESSEL(X)/X
double fncBx(const double x)
{
    double valSign= 1.;
     double sum = 0.5;
     double y =  x* x /4.;
     double prod = 0.5;
     for (int i = 1; i < 50; ++i)
     {
        prod = -prod *  y/ ((double) (i * (i +1))) ;
        sum += prod;
     }


     return sum;
}

// производная фукнкции BESSEL(X)/X
double fnc_dBx_po_dx(const double x)
{
    double valSign= -1.;
    double sum = - x/ 4.;
    double y = x* x/4.;
    double prod = - x/ 4.;
    for (int i = 2; i < 40; ++i)
    {

       prod = -prod *  y/ ((double) (i * (i +1))) ;
       sum += prod * ((double) i);
    }
    return sum/2.;
}
//----------------------------------
double QReceiver::calcSigmaNoise()
{
    return 0.00024/sqrt(2.);
}
//-------------------------------------------
void QReceiver::collectData_forDiagrGraph(double lamb, int QUantRows
                                          , double step, double *arrBuff, double *arrBuffPolar)
{

int QUantColsReport0 = 2;
for (int i =0; i < QUantRows; ++i)
{

    double tetta = -M_PI/2. + ((double)i) * step;
    arrBuff[i *QUantColsReport0 ] = tetta ;
    arrBuff[i *QUantColsReport0 +1] = fncBessDiagr(tetta,lamb);

    arrBuffPolar[i *QUantColsReport0 ] = arrBuff[i *QUantColsReport0 +1]* sin(tetta);
    arrBuffPolar[i *QUantColsReport0 +1] = arrBuff[i *QUantColsReport0 +1]* cos(tetta);

}
}
