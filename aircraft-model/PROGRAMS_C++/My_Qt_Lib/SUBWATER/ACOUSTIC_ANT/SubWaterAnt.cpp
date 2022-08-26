#include "SubWaterAnt.h"
#include <math.h>
#include "Comp.h"
#include "Receiver.h"
#include "Gauss.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

using namespace std;

QSubWaterAnt::QSubWaterAnt():
    mQuant(0)
   ,ma(0)

{

}
// конструктор копирования
 QSubWaterAnt ::QSubWaterAnt (const QSubWaterAnt &R)
 {
     mvctReceiver = R.mvctReceiver;
     mvctPntXY = R.mvctPntXY;

     mQuant = R.mQuant;
     ma = R.ma;

 }


// оператор присваивания
QSubWaterAnt &QSubWaterAnt::operator=(const QSubWaterAnt  &R)
{
    mvctReceiver = R.mvctReceiver;
    mvctPntXY = R.mvctPntXY;

    mQuant = R.mQuant;
    ma = R.ma;

  return *this ;
}
// парам констр
QSubWaterAnt :: QSubWaterAnt( const QVector<QReceiver> vctSubWaterAnt, const QVector<QPointF> vctPntXY)

{
   mvctReceiver = vctSubWaterAnt ;
   mvctPntXY = vctPntXY;
}

//----------------------------------------
// парам констр
 QSubWaterAnt :: QSubWaterAnt(const int QuantReceivers, const double a_Receiver
                                  , const double radius)

{
   mvctPntXY.resize(QuantReceivers);
   mvctReceiver.resize(QuantReceivers);

   QReceiver Receiver (a_Receiver);
   mvctReceiver.fill(Receiver,QuantReceivers);
   mQuant = QuantReceivers;

  ma = a_Receiver;

}


//---------------------------------------
TComp QSubWaterAnt ::calcDiagrB(const double lamb,const double q, const double e)
{
 // единичный вектор нормали к плоскости волны
    double arr_e[3]= {0.};
    arr_e[0] = -cos (e) * sin(q);
    arr_e[1] = -cos (e) * cos(q);
    arr_e[2] = -sin (e) ;
    // !

    // массив набега фаз
    // вычисляется расстояние от точки mvctPntXY[i] до плоскости,
    // которая задается тчкой (0,0,0) и вектором нормали arr_e
    double valNabeg = 0.;
    double transmitterWaveAmp = 1.;
    double dist = 1.;

    double valk = 2. * M_PI/ lamb;

     TComp cmpDiagr(0.,0.);
     for (int i =0; i < mvctPntXY.size(); ++i)
    {
      valNabeg = valk*(mvctPntXY.at(i).x() * arr_e[0] + mvctPntXY.at(i).y() * arr_e[1]);

      TComp cmp0 (cos(valNabeg), sin (valNabeg));

      QReceiver receiver = mvctReceiver.at(i);
      double amp = receiver.fncBessDiagr(M_PI /2. - abs(e), lamb);

      cmp0 *=  TComp(amp, 0) ;
      cmpDiagr = cmpDiagr + cmp0;

    }
     cmpDiagr *= TComp(1./ ((double)mvctPntXY.size()), 0.);
     return cmpDiagr;
}
//------------------------------------

double QSubWaterAnt ::calcF(double c1, double c2, double lam1, double lam2, double  al0)
{
    return 1. - c1 * c1 /((lam1+al0) *(lam1+al0))- c2 * c2 /((lam2+al0) *(lam2+al0));
}

double QSubWaterAnt ::calcdF(double c1, double c2, double lam1, double lam2, double  al0)
{
    return 2.*c1 * c1 /((lam1+al0) *(lam1+al0)*(lam1+al0)) +2.* c2 * c2 /((lam2+al0) *(lam2+al0)*(lam2+al0));
}
//----------------------------------
void QSubWaterAnt::calcVectMeasures(const double lamb,double valSignalAmp, double q, double e, TComp *arrQ )
{
    // единичный вектор нормали к плоскости волны
       double arr_e[3]= {0.};
       arr_e[0] = -cos (e) * sin(q);
       arr_e[1] = -cos (e) * cos(q);
       arr_e[2] = -sin (e) ;
       // !

       // массив набега фаз
       // вычисляется расстояние от точки mvctPntXY[i] до плоскости,
       // которая задается тчкой (0,0,0) и вектором нормали arr_e

      // double valNabeg = 0.;
       double valk = 2. * M_PI/ lamb;

        TComp cmpDiagr(0.,0.);
        double *arrNabegL = new double [mvctPntXY.size()];
        double *arrNabegPh= new double [mvctPntXY.size()];
        for (int i =0; i < mvctPntXY.size(); ++i)
       {
            QReceiver Receiver = mvctReceiver.at(i);
         double diagr =  Receiver.fncBessDiagr( M_PI/2. - fabs(e), lamb);
         arrNabegL[i] = (mvctPntXY.at(i).x() * arr_e[0] + mvctPntXY.at(i).y() * arr_e[1]);
         arrNabegPh[i] = valk*arrNabegL[i];
         TComp cmp0 (valSignalAmp * diagr*cos(arrNabegPh[i]),valSignalAmp * diagr* sin(arrNabegPh[i]));

         arrQ [i]= cmp0;
       }

}

//
// q0, tetta_scan - куда ориентрована диаграмма
// по углу места относительо нормали к антенне
// e - угол места относительноплоскости антенны
TComp QSubWaterAnt ::calcSumNormalizedDiagr(const double lamb,const double q0, const double tetta_scan
                                ,const double q, const double e)
{
 // единичный вектор нормали направления q0,e0
    double arr_e[3]= {0.};
    arr_e[0] = -cos (e) * sin(q);
    arr_e[1] = -cos (e) * cos(q);
    arr_e[2] = -sin (e) ;
    // !

    // массив набега фаз
    // вычисляется расстояние от точки mvctPntXY[i] до плоскости,
    // которая задается тчкой (0,0,0) и вектором нормали arr_e
    double valNabeg = 0.;   
    double valk = 2. * M_PI/ lamb;
     TComp cmpDiagr(0.,0.);
     TComp *cmparrMultipliers = new TComp[mvctPntXY.size()];
     calcMultiplierArr( lamb,  q0,  M_PI/2. -tetta_scan,cmparrMultipliers);

     for (int i =0; i < mvctPntXY.size(); ++i)
    {
     QReceiver receiver = mvctReceiver.at(i);
     double amp = receiver.fncBessDiagr(M_PI /2. - abs(e), lamb);

      valNabeg = valk *(mvctPntXY.at(i).x() * arr_e[0] + mvctPntXY.at(i).y() * arr_e[1]);

      TComp cmp0 (cos(valNabeg), sin (valNabeg));
      cmp0 *=  cmparrMultipliers[i] ;

      cmp0 *= TComp( amp, 0.);
      cmpDiagr = cmpDiagr + cmp0;

    }
     cmpDiagr *= TComp(1./ ((double)mvctPntXY.size()), 0.);
     delete []cmparrMultipliers;
     return cmpDiagr;
}
//----------------------------------
// вычисление множителей формирования диаграммы с направления q0, e0
void QSubWaterAnt ::calcMultiplierArr(const double lamb, const double q0, const double e0
                                   ,TComp *cmparrMultipliers)
{
    // единичный вектор нормали направления q0,e0
       double arr_e0[3]= {0.};
       arr_e0[0] = -cos (e0) * sin(q0);
       arr_e0[1] = -cos (e0) * cos(q0);
       arr_e0[2] = -sin (e0) ;
       // !
       double valk = 2. * M_PI/ lamb;
       for (int i =0; i < mvctPntXY.size(); ++i)
      {
        double valPhNabeg = valk*(mvctPntXY.at(i).x() * arr_e0[0]
                + mvctPntXY.at(i).y() * arr_e0[1]);
        cmparrMultipliers[i] = TComp(cos(valPhNabeg), -sin(valPhNabeg));

      }
}
//--------------------------
//------------------------------------
void QSubWaterAnt ::calcAmpRaznMeth_e(const double lamb,TComp *arrQ, double q0, double e0, double *pe)
{
    // вычисление ширины диагнраммы влево и вправо и производныз
    double tetta_scan0 = M_PI/2. - e0;
    TComp cmp_amp0 = calcSumNormalizedDiagr( lamb, q0,  tetta_scan0, q0,  e0);
    double amp0= cmp_amp0.modul();
    double diagr_prev = amp0;
    double step  = 0.0001;
    double d_diagr_r = 0.;
    double tet_07 = 0.;
    for (int i= 0; i < 10000; ++i)
    {
        tet_07 = ((double) i) * step;
        TComp cmp_amp = calcSumNormalizedDiagr( lamb, q0,  tetta_scan0, q0,  e0 + tet_07);
        double diagr= cmp_amp.modul();        

        if (diagr <= 0.7 * amp0)
        {
            break;
        }

    }


    // !
    double tet_0 = 0;


    double der1 = 0., der2 =2.;

    double diagr1 =0., diagr2 = 0.;
    for (int i=0; i< 100000; ++i)
    {
       tet_0 = e0 - tet_07 +((double) i) * step;
       diagr1 =  calcSumNormalizedDiagr( lamb, q0, M_PI/2. -(e0 - tet_07), q0,  tet_0).modul();
       diagr2 =  calcSumNormalizedDiagr( lamb, q0, M_PI/2. -( e0 + tet_07), q0,  tet_0).modul();
       if (fabs(diagr2 -diagr1)< 0.001)
       {
          //der1 = (t1 - t1_pr) /step;
         // der2 = (t2 - t2_pr) /step;

          double x1= calcSumNormalizedDiagr( lamb, q0,  e0 - tet_07, q0,  tet_0+0.001).modul();
          double x2 =calcSumNormalizedDiagr( lamb, q0,  e0 - tet_07, q0,  tet_0).modul();
                  double x3 = calcSumNormalizedDiagr( lamb, q0,  e0 + tet_07, q0,  tet_0 + 0.001).modul();
                  double x4 =calcSumNormalizedDiagr( lamb, q0,  e0 + tet_07, q0,  tet_0).modul();
          der1= (calcSumNormalizedDiagr( lamb, q0,  M_PI/2. -(e0 - tet_07), q0,  tet_0+0.001).modul()
                                 - calcSumNormalizedDiagr( lamb, q0,  M_PI/2. -(e0 - tet_07), q0,  tet_0).modul())/0.001;

          der2= (calcSumNormalizedDiagr( lamb, q0,  M_PI/2. -( e0 + tet_07), q0,  tet_0 + 0.001).modul()
                          -calcSumNormalizedDiagr( lamb, q0, M_PI/2. -( e0 + tet_07), q0,  tet_0).modul())/0.001;
          break;
       }

    }

    ////


    TComp cmp_diagr_right = calcIzm(  lamb,arrQ, q0,  e0 -tet_07);
    double aright = cmp_diagr_right.modul() ;

    TComp cmp_diagr_left = calcIzm(  lamb,arrQ, q0,  e0 +tet_07);
    double aleft = cmp_diagr_left.modul() ;
    double ro = (aright - aleft)/(aright + aleft);
    double mu = ro * (der1 + der2);

    double del_e = (2. *diagr1 * ro)/ (der1 -  der2 - mu);



    *pe = tet_0 + del_e;


}

//------------------------------------
void QSubWaterAnt ::calcAmpRaznMeth_q(const double lamb,TComp *arrQ, double q0, double e0, double *pq)
{
    // вычисление ширины диагнраммы влево и вправо и производныз
    double tetta_scan0 = M_PI/2. - e0;
    TComp cmp_amp0 = calcSumNormalizedDiagr( lamb, q0,  tetta_scan0, q0,  e0);
    double amp0= cmp_amp0.modul();
    double diagr_prev = amp0;
    double step  = 0.0001;
    double d_diagr_r = 0.;
    double tet_07 = 0.;
    for (int i= 0; i < 10000; ++i)
    {
        tet_07 = ((double) i) * step;
        TComp cmp_amp = calcSumNormalizedDiagr( lamb, q0,  tetta_scan0, q0+ tet_07,  e0 );
        double diagr= cmp_amp.modul();

        if (diagr <= 0.7 * amp0)
        {
            break;
        }

    }


    // !
    double q_0 = 0;


    double der1 = 0., der2 =2.;

    double diagr1 =0., diagr2 = 0.;
    for (int i=0; i< 100000; ++i)
    {
       q_0 = q0 - tet_07 +((double) i) * step;
       diagr1 =  calcSumNormalizedDiagr( lamb, q0 +tet_07, tetta_scan0, q_0,  e0).modul();
       diagr2 =  calcSumNormalizedDiagr( lamb, q0-tet_07, tetta_scan0, q_0,  e0).modul();
       if (fabs(diagr2 -diagr1)< 0.001)
       {

          der1= (calcSumNormalizedDiagr( lamb, q0 +tet_07, tetta_scan0, q_0 + 0.001,  e0).modul()
                                 - calcSumNormalizedDiagr( lamb, q0 +tet_07, tetta_scan0, q_0,  e0).modul())/0.001;

          der2= (calcSumNormalizedDiagr( lamb, q0 -tet_07, tetta_scan0, q_0 + 0.001,  e0).modul()
                 - calcSumNormalizedDiagr( lamb, q0 -tet_07, tetta_scan0, q_0,  e0).modul())/0.001;
          break;
       }

    }

    ////


    TComp cmp_diagr_right = calcIzm(  lamb,arrQ, q0-tet_07,  e0 );
    double aright = cmp_diagr_right.modul() ;

    TComp cmp_diagr_left = calcIzm(  lamb,arrQ, q0+tet_07,  e0 );
    double aleft = cmp_diagr_left.modul() ;
    double ro = (aright - aleft)/(aright + aleft);
    double mu = ro * (der1 + der2);

    double del_e = (2. *diagr1 * ro)/ (der1 -  der2 - mu);



    *pq = q_0 + del_e;


}


TComp QSubWaterAnt ::calcIzm( const double lamb,TComp *arrQ,const double q0, const double e0 )
{
    TComp *cmparrMultipliers = new TComp[mvctPntXY.size()];
    calcMultiplierArr(lamb,  q0,  e0,cmparrMultipliers);
    TComp rez(0.,0.);
    for(int i =0; i < mvctPntXY.size(); ++i)
    {
       TComp temp =  arrQ[i] * cmparrMultipliers[i];
       rez += temp;
    }

    delete[]cmparrMultipliers;
    return rez;
}
//--------------------------------
// набивается массив arrBuff[QUantCols * QUantRows]
//q0, e0 - направление куда смотрит диаграмма
//step - шаг по углу e с которым вычисляется график
// выводится амплитуда и фаза (аргумент)
void QSubWaterAnt ::collectData_for_e_DiagrGraph(const double lamb,const double q0
                  ,const double e0, const int QUantRows,const double step,double *arrBuff)
{
    const int QUantCols =3;
    for (int i =0; i < QUantRows; ++i)
         {
             double e = step + ((double)i) * step;

             arrBuff[i *QUantCols ] =e;

             TComp Comp1 = calcSumNormalizedDiagr(lamb,q0, e0
                                             ,q0,e);

                     arrBuff[i *QUantCols +1] = Comp1.modul();
                     arrBuff[i *QUantCols +2] = Comp1.phase();
         }
}
//---------------------------
void QSubWaterAnt ::collectData_forDiagrGraph_e(double lamb,double q0,double tetta0 , int QUantRows
                                    ,double  step,double * arrBuff,double * arrBuffPolar)
{

    int QUantCols = 3;
for (int i =0; i < QUantRows; ++i)
{

    double e =  step + ((double)i) * step;

    arrBuff[i *QUantCols ] = M_PI/2. - fabs(e);

    TComp Comp = calcSumNormalizedDiagr(lamb,q0, tetta0
                                    ,q0,e);
            arrBuff[i *QUantCols +1] = Comp.modul();
            arrBuff[i *QUantCols +2] = Comp.phase();

            arrBuffPolar[i *2] = arrBuff[i *QUantCols +1] * cos(e);
            arrBuffPolar[i *2 +1] = arrBuff[i *QUantCols +1] * sin(e);

}
}

//--------------------------------
void QSubWaterAnt ::fillReceiversLocation(const double radius)
{
int ii=0;
}


