#include "RingedAnt.h"
#include <math.h>
#include "Comp.h"
#include "Receiver.h"
#include "Gauss.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

using namespace std;

QRingedAnt::QRingedAnt():
     mr(0)
    ,mQuant(0)
    ,ma(0)
    ,mscatter_a(0)
    ,mscatter_ang(0)
    ,mscatter_radius(0)
{

}
// конструктор копирования
 QRingedAnt ::QRingedAnt (const QRingedAnt &R)
 {
     mvctReceiver = R.mvctReceiver;
     mvctPntXY = R.mvctPntXY;
     mr = R.mr;
     mQuant = R.mQuant;
     ma = R.ma;
     mscatter_a = R.mscatter_a;
     mscatter_ang = R.mscatter_ang;
     mscatter_radius  = R.mscatter_radius ;
 }


// оператор присваивания
QRingedAnt &QRingedAnt::operator=(const QRingedAnt  &R)
{
    mvctReceiver = R.mvctReceiver;
    mvctPntXY = R.mvctPntXY;
    mr = R.mr;
    mQuant = R.mQuant;
    ma = R.ma;
    mscatter_a = R.mscatter_a;
    mscatter_ang = R.mscatter_ang;
    mscatter_radius  = R.mscatter_radius ;
  return *this ;
}
// парам констр
QRingedAnt :: QRingedAnt( const QVector<QReceiver> vctRingedAnt, const QVector<QPointF> vctPntXY)

{
   mvctReceiver = vctRingedAnt ;
   mvctPntXY = vctPntXY;
}

//----------------------------------------
// парам констр
 QRingedAnt :: QRingedAnt(const int QuantReceivers, const double a_Receiver
                                  , const double radius)

{
   mvctPntXY.resize(QuantReceivers);
   mvctReceiver.resize(QuantReceivers);

   QReceiver Receiver (a_Receiver);
   mvctReceiver.fill(Receiver,QuantReceivers);


   double step = 2. * M_PI/ ((double)QuantReceivers);

   for (int i = 0; i < QuantReceivers; ++i)
   {
       double ang = ((double)i) * step;
       mvctPntXY.replace(i, QPointF(radius * sin(ang),radius * cos(ang)));
   }
  mQuant = QuantReceivers;
  ma = a_Receiver;
  mr = radius;

}

 // парам констр
  QRingedAnt :: QRingedAnt(const int QuantReceivers, const double a
            , const double radius, const double scatter_a
                         , const double scatter_ang, const double scatter_radius )

 {
    mvctPntXY.resize(QuantReceivers);
    mvctReceiver.resize(QuantReceivers);

    for (int i =0; i < QuantReceivers; ++i)
    {
       QReceiver Receiver (a -scatter_a/ 2. + scatter_a * getRand01());
       mvctReceiver.replace(i,Receiver );
    }

    double step = 2. * M_PI/ ((double)QuantReceivers);

    for (int i = 0; i < QuantReceivers; ++i)
    {
        double ang =  ((double)i) * step - scatter_ang/2 + scatter_ang * getRand01();
        double r = radius - scatter_radius/2. + scatter_radius * getRand01();
        mvctPntXY.replace(i, QPointF(r * sin(ang),r * cos(ang)));
    }

    mr = radius;
    mQuant = QuantReceivers;
    ma = ma;
    mscatter_a = scatter_a;
    mscatter_ang = scatter_ang;
    mscatter_radius  = scatter_radius ;
 }
//---------------------------------------
TComp QRingedAnt ::calcDiagrB(const double lamb,const double q, const double e)
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
void QRingedAnt ::calcAngs(const double lamb,TComp *arrQ, double *pq, double *pe)
{
    // формирование матьрицы a
    double *pmtrx_a = new double[(mQuant - 1) *2];
    memset(pmtrx_a, 0, sizeof(double) *(mQuant - 1) *2);
    for (int i =0; i < ((mQuant - 1)); ++i)
    {
      pmtrx_a[2 * i] =  mvctPntXY.at(i + 1).x() - mvctPntXY.at(0).x()  ;
      pmtrx_a[2 * i + 1] =  mvctPntXY.at(i + 1).y() - mvctPntXY.at(0).y() ;

     // double t0 = mr * sin(2. *M_PI/ ((double)mQuant) * ((double)(i + 1.)));
     //.
     // double t1 = -mr * 2. *sin(M_PI/ ((double)mQuant) * ((double)(i + 1.)))*sin(M_PI/ ((double)mQuant) * ((double)(i + 1.)));

     //                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  int yy=0;
    }
    // !

    // формирование вектора b
    double fi0 = arrQ[0].phase();
    double *parr_b= new double[(mQuant - 1)];
    memset(parr_b, 0, sizeof(double) *(mQuant - 1) );
    double valk = lamb/2./M_PI;
    for (int i =0; i < ((mQuant - 1)); ++i)
    {
      parr_b[i] = valk *(arrQ[i + 1].phase() -fi0);
    }
    // !

    // вычисляем вектор с[2]
    double arrc[2] = {0.};
    MtrxTranspMultMatrx( pmtrx_a,mQuant - 1, 2, parr_b,1, arrc) ;
    // !

    // решение ураынения
    double val_quant =(double) mQuant;
    double lam1 = mr * mr * (val_quant /2.);
    double lam2 = 4. *( (val_quant -1.) * 3.
                        + 4. * tan(M_PI/ val_quant) - tan(2 *M_PI/ val_quant));

    double arr_ata[4] ={0.}, arrTemp[100];
    memcpy(arrTemp,pmtrx_a, 2 * (mQuant - 1) * sizeof(double));
    MtrxTranspMultMatrx(pmtrx_a,mQuant - 1, 2, arrTemp,2, arr_ata) ;
    double al0 =(lam1 + lam2)/ 2.;

    lam2 = arr_ata[3];
    lam1 = arr_ata[0];
    al0 = 0.00001;

    double temp = arrc[0] * arrc[0] /lam1 /lam1 + arrc[1] * arrc[1] /lam2 /lam2;

    for ( int i =0; i <200; ++i)
    {
        double f = calcF(arrc[0], arrc[1], lam1, lam2, al0);
        double dF = calcdF(arrc[0], arrc[1], lam1, lam2, al0);
        double del =f/ dF;
        al0 -= 0.5 * del;
        if (fabs(del)< 0.00001)
        {
            break;
        }
    }
    // !

    double x1 = arrc[0]/ (arrc[0] + al0);
    double x2 = arrc[1]/ (arrc[1] + al0);

   *pe =asin(1.-x1 * x1);
    *pq =atan2(x1,x2);
    delete []pmtrx_a ;
    delete []parr_b;


    // единичный вектор нормали к плоскости волны
    // ОТЛАДКА
    double arr_e[3]= {0.};
    double e = M_PI/ 4., q = M_PI/4.;
    arr_e[0] = cos (e) * sin(q);
    arr_e[1] = cos (e) * cos(q);
    arr_e[2] = sin (e) ;

    double al1 = arrc[0]/ arr_e[0]  -lam1;
    double al2 = arrc[1]/ arr_e[1]  -lam2;
    int rr = 0;
    // ОТЛАДКА !


}
double QRingedAnt ::calcF(double c1, double c2, double lam1, double lam2, double  al0)
{
    return 1. - c1 * c1 /((lam1+al0) *(lam1+al0))- c2 * c2 /((lam2+al0) *(lam2+al0));
}

double QRingedAnt ::calcdF(double c1, double c2, double lam1, double lam2, double  al0)
{
    return 2.*c1 * c1 /((lam1+al0) *(lam1+al0)*(lam1+al0)) +2.* c2 * c2 /((lam2+al0) *(lam2+al0)*(lam2+al0));
}
//----------------------------------
void QRingedAnt::calcVectMeasures(const double lamb,double valSignalAmp, double q, double e, TComp *arrQ )
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
TComp QRingedAnt ::calcSumNormalizedDiagr(const double lamb,const double q0, const double tetta_scan
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
void QRingedAnt ::calcMultiplierArr(const double lamb, const double q0, const double e0
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
void QRingedAnt ::calcAmpRaznMeth_e(const double lamb,TComp *arrQ, double q0, double e0, double *pe)
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
void QRingedAnt ::calcAmpRaznMeth_q(const double lamb,TComp *arrQ, double q0, double e0, double *pq)
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


TComp QRingedAnt ::calcIzm( const double lamb,TComp *arrQ,const double q0, const double e0 )
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
void QRingedAnt ::collectData_for_e_DiagrGraph(const double lamb,const double q0
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
void QRingedAnt ::collectData_forDiagrGraph_e(double lamb,double q0,double tetta0 , int QUantRows
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


