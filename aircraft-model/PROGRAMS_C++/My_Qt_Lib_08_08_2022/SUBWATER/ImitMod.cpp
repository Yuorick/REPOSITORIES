#include "ImitMod.h"
#include "MatrixProccess.h"
#include <math.h>
#include <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "Gauss.h"
#include "SubWaterBeam.h"
#include "TrueMeasParams.h"
#include "CoordSystTrsf.h"


#include "BigMeasure.h"


class QBigMeasure;
QImitMod::QImitMod()
{
     mVess = QPeaceVess();

     mtblRealPrfl = TTable_1D();

     memset(marrTrueBeaconPos, 0, 3 * sizeof(double));
     memset(marrAprioriBeaconPos, 0, 3 * sizeof(double));
     memset(marrTrueBeaconVelo, 0, 3 * sizeof(double));
     memset(marrAprioriBeaconVelo, 0, 3 * sizeof(double));

     mTypeOfVessTraj = PNT_6;
     mTObrabotki = 0.;
     memset(marrAprioriGpsPosParams, 0, 6 * sizeof(double));


}

// конструктор копирования
  QImitMod :: QImitMod (const  QImitMod &R)
 {
   mVess = R.mVess;
   mtblRealPrfl = R.mtblRealPrfl;
   memcpy(marrTrueBeaconPos, R.marrTrueBeaconPos, 3 * sizeof(double));
   memcpy(marrAprioriBeaconPos, R.marrAprioriBeaconPos, 3 * sizeof(double));
   memcpy(marrTrueBeaconVelo, R.marrTrueBeaconVelo, 3 * sizeof(double));
   memcpy(marrAprioriBeaconVelo, R.marrAprioriBeaconVelo, 3 * sizeof(double));
   mTypeOfVessTraj = R.mTypeOfVessTraj;
   mTObrabotki = R.mTObrabotki;
   memcpy(marrAprioriGpsPosParams, R.marrAprioriGpsPosParams, 6 * sizeof(double));
 }

 // оператор присваивания
 QImitMod  &QImitMod::operator=( const QImitMod  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      mVess = R.mVess;
      mtblRealPrfl = R.mtblRealPrfl;
      memcpy(marrTrueBeaconPos, R.marrTrueBeaconPos, 3 * sizeof(double));
      memcpy(marrAprioriBeaconPos, R.marrAprioriBeaconPos, 3 * sizeof(double));
      memcpy(marrTrueBeaconVelo, R.marrTrueBeaconVelo, 3 * sizeof(double));
      memcpy(marrAprioriBeaconVelo, R.marrAprioriBeaconVelo, 3 * sizeof(double));
      mTypeOfVessTraj = R.mTypeOfVessTraj;
      mTObrabotki = R.mTObrabotki;
      memcpy(marrAprioriGpsPosParams, R.marrAprioriGpsPosParams, 6 * sizeof(double));
     return *this ;
 }

 // парам конструктор 1
QImitMod:: QImitMod (const QPeaceVess Vess, const TTable_1D tblRealPrfl
           ,const double *arrTrueBeaconPos,const double *arrAprioriBeaconPos
           , const enumTypeOfVessTraj TypeOfVessTraj, const double TObrabotki
                     ,const double *arrAprioriGpsPosParams)

{
   mVess = Vess;
   mtblRealPrfl = tblRealPrfl;
   mTypeOfVessTraj = TypeOfVessTraj;
   memcpy(marrTrueBeaconPos, arrTrueBeaconPos, 3 * sizeof(double));
   memcpy(marrAprioriBeaconPos, arrAprioriBeaconPos, 3 * sizeof(double));
   memset(marrTrueBeaconVelo, 0, 3 * sizeof(double));
   memset(marrAprioriBeaconVelo, 0, 3 * sizeof(double));
   mTObrabotki = TObrabotki;
   memcpy(marrAprioriGpsPosParams, arrAprioriGpsPosParams, 6 * sizeof(double));
}
//---------------------------------
// парам конструктор 2
QImitMod:: QImitMod (const QPeaceVess Vess, const TTable_1D tblRealPrfl
          ,const double *arrTrueBeaconPos,const double *arrAprioriBeaconPos
          , const enumTypeOfVessTraj TypeOfVessTraj, const double TObrabotki
       ,const double *arrAprioriGpsPosParams, const double *arrTrueBeaconVelo, const double *arrAprioriBeaconVelo )

{
  mVess = Vess;
  mtblRealPrfl = tblRealPrfl;
  mTypeOfVessTraj = TypeOfVessTraj;
  memcpy(marrTrueBeaconPos, arrTrueBeaconPos, 3 * sizeof(double));
  memcpy(marrAprioriBeaconPos, arrAprioriBeaconPos, 3 * sizeof(double));
  memcpy(marrTrueBeaconVelo, arrTrueBeaconVelo, 3 * sizeof(double));
  memcpy(marrAprioriBeaconVelo, arrAprioriBeaconVelo, 3 * sizeof(double));
  mTObrabotki = TObrabotki;
  memcpy(marrAprioriGpsPosParams, arrAprioriGpsPosParams, 6 * sizeof(double));
}
 //-------------------------------------------------
 void QImitMod::createMeasures ( const int QuantMeas, QBigMeasure *arrMeas)
 {
     const double DeepthMax = mtblRealPrfl.mparrArg[mtblRealPrfl.mNumCols - 2];

     double ZonaR_ = calcXMCrit( 2., DeepthMax, mtblRealPrfl);
     if (ZonaR_ > 1000.)
     {
         ZonaR_ = 1000.;
     }
     const double VAlQdrtSideLength = ZonaR_* 0.4 *2.;//ZonaR_* 0.66 *2.;

     switch(mTypeOfVessTraj)
     {
     case PNT_6:
     case PNT_7:
         createSeparatedPointsMeasures(QuantMeas,ZonaR_,arrMeas);
         break;

     case QDRT:
         createQdrtTrajMeasures(QuantMeas,ZonaR_,arrMeas);
         break;

     case ZIG_ZAG:
       createZigZagTrajMeasures(QuantMeas,ZonaR_,arrMeas);
        break;

     case DIAMETRS:
       createDiamTrajMeasures(QuantMeas,ZonaR_,arrMeas);
        break;

     case LINE6:
       createLine6TrajMeasures(QuantMeas,ZonaR_,arrMeas);
        break;

     default:
         break;
     }

     /*if (mTypeOfVessTraj < 2)
     {
        for (int i =0;i < QuantMeas-1;++i)
        {
        // bear - пеленг корабля относительно маяка
        // double bear = -2. * M_PI/ (2.*((double)QuantMeas)) + 2. * M_PI/((double)QuantMeas) *((double)i);
        double bear =  2. * M_PI/((double)QuantMeas) *((double)i);
        // курсовой угол скорости корабля
        mVess.mQ0 = bear + M_PI/2.;

        mVess.marrVectSost[0] = marrAprioriBeaconPos[0] +  sin(bear) *ZonaR_ *0.4;
        mVess.marrVectSost[1] = marrAprioriBeaconPos[1] +  cos(bear) *ZonaR_ *0.4;
        mVess.marrVectSost[3]= 0.;
        mVess.marrVectSost[4]= 0.;
        mVess.marrVectSost[5]= 0.;
        mVess.mVQ = 0.;
        mVess.mVVess = 0.;
        mVess.mQ = mVess.mQ0; // угол курса
        mVess.mPsi = 0.; // угол килевой качки
        mVess.mTet = 0.; // угол бортовой качки
        mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];

        //memcpy(mVess.marrVectSost0, mVess.marrVectSost,9 * sizeof(double));
        mVess.tuneCurrentSinsMeasures();
        createSingleMeasure(arrMeas[i]);

        }

        // по центру

        mVess.mQ0 =  M_PI/2.;

        mVess.marrVectSost[0] = marrAprioriBeaconPos[0] ;
        mVess.marrVectSost[1] = marrAprioriBeaconPos[1] ;
        mVess.marrVectSost[3]= 0.;
        mVess.marrVectSost[4]= 0.;
        mVess.marrVectSost[5]= 0.;
        mVess.mVQ = 0.;
        mVess.mVVess = 0.;
        mVess.mQ = mVess.mQ0; // угол курса
        mVess.mPsi = 0.; // угол килевой качки
        mVess.mTet = 0.; // угол бортовой качки
        mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];

        //memcpy(mVess.marrVectSost0, mVess.marrVectSost,9 * sizeof(double));
        mVess.tuneCurrentSinsMeasures();
        createSingleMeasure(arrMeas[QuantMeas-1]);
     }
     else
     {

         int numQdrtSide = 0;
         initPartOfQdrtTraj(numQdrtSide,marrAprioriBeaconPos, VAlQdrtSideLength);

         for (int i =0; i < QuantMeas; ++i)
         {
           int numQdrtSideCur = calcNumQdrtSide(marrAprioriBeaconPos,  VAlQdrtSideLength
                                                , numQdrtSide);
           if(numQdrtSideCur != numQdrtSide)
           {
             numQdrtSide = numQdrtSideCur ;
             initPartOfQdrtTraj(numQdrtSide,marrAprioriBeaconPos, VAlQdrtSideLength);
           }

           //pln.Points[i].X =mVess.marrVectSost[0];
           //pln.Points[i].Y =mVess.marrVectSost[1];
           createSingleMeasure(arrMeas[i]);
         }


     }
*/

 }

//--------------------------------------------
int QImitMod::calcNumQdrtSide(const double *arrQdrtCentre, const double VAlQdrtSideLength
                              , const int NUmQdrtSideCur)
{
    // координаты корабля в плоскости XY
    double valPntX = mVess.marrVectSost[0];
    double valPntY = mVess.marrVectSost[1];

    // координаты квадраты по оси OX
    double ax = arrQdrtCentre[0] - VAlQdrtSideLength/2;
    double bx = arrQdrtCentre[0] + VAlQdrtSideLength/2;

    // координаты квадраты по оси OY
    double ay = arrQdrtCentre[1] - VAlQdrtSideLength/2;
    double by = arrQdrtCentre[1] + VAlQdrtSideLength/2;

    // условие того, что корабль находится на стороне квадрата

    if( (fabs(valPntX - arrQdrtCentre[0]) <= (VAlQdrtSideLength/2. + 0.1) )
     &&(fabs(valPntY - arrQdrtCentre[1]) <= (VAlQdrtSideLength/2. + 0.1) ))
    {
        return NUmQdrtSideCur;
    }
    return (NUmQdrtSideCur + 1) % 4 ;
}
 //-------------------------------------------------
 void  QImitMod::initPartOfQdrtTraj(const int numQdrtSide
                 ,const double *arrAprioriBeaconPos,const double VAlQdrtSideLength)
 {
 // bear - пеленг корабля относительно маяка
 double bear = - M_PI/ 4. +  M_PI/2. *((double)numQdrtSide);
 // курсовой угол скорости корабля
 mVess.mQ0 = M_PI/2. + M_PI/2. *((double)numQdrtSide);
 mVess.mQ = mVess.mQ0; // угол курса

 double hipoten = VAlQdrtSideLength/ 2.*sqrt(2.);

 mVess.marrVectSost[0] = marrAprioriBeaconPos[0] +  sin(bear) *hipoten;
 mVess.marrVectSost[1] = marrAprioriBeaconPos[1] +  cos(bear) *hipoten;
 mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];
 mVess.marrVectSost[3]= mVess.mVVess * sin(mVess.mQ0);
 mVess.marrVectSost[4]= mVess.mVVess * cos(mVess.mQ0);
 mVess.marrVectSost[5]= 0.;
 mVess.mVQ = 0.;

 mVess.tuneCurrentSinsMeasures();
 }

//-----------------------------------
 bool QImitMod::createSingleMeasure ( QBigMeasure &BigMeas)
 {
    QTrueMeasParams trueMeasParams;
    if(!mVess.createTrueMeasureParams( mtblRealPrfl, mTObrabotki
                                   , marrTrueBeaconPos,trueMeasParams))
    {
        return false;
    }
    imitateMeasure(trueMeasParams, BigMeas);

    // ОТЛАДКА
   // double arrSgsk[3]={0.}, valTZv = -1;;
    //transfMeasure_to_GSK(mtblRealPrfl,BigMeas,mVess.mPlatform.marrPosParams,/*marrTrueBeacVelo*/NULL,  arrSgsk,   &valTZv);

    // !ОТЛАДКА
    return true;
 }
 //------------------------------------------
 // имитация измерения антенны
 // INPUT:
 //trueMeasParams - структура истинных параметров, харакетеризующих положение цели
 // и обеспечивающих возможность наложения шумов и погрешностей
 //OUTPUT:
 //Meas - структура, в которой хранятся данные одного измерения
void QImitMod:: imitateMeasure(QTrueMeasParams &trueMeasParams, QBigMeasure &BigMeas)
{
    // 1. формирование измерений антенны
    mVess.mPlatform.mpHidroRLS->imitateMesure( trueMeasParams, BigMeas);
   /* double valTzaprZv,  valTotvZv  = 0.,  valSig_t  = 0.,  val_qzv
                                         = 0.,  val_ezv  = 0.,  valSig_q  = 0.,  valSig_e =0.;
    mVess.mPlatform.mpHidroRLS->imitateMesure( trueMeasParams,    &valTzaprZv
                                   ,    &valTotvZv,    &valSig_t,    &val_qzv
                                   ,    &val_ezv,    &valSig_q,    &valSig_e);*/

    // 1!

    // 2. формирование измерения положения судна на момент запросного сигнала
      imitateVesselPosition(trueMeasParams.marrSVessWave,trueMeasParams.marrMuWave
                    ,trueMeasParams.marrMuWaveZv,BigMeas.marrSVessWaveZv);
      // 2!

      // 3. формирование измерения положения судна на момент ответного сигнала
        imitateVesselPosition(trueMeasParams.marrSVess,trueMeasParams.marrMu
                      ,trueMeasParams.marrMuZv,BigMeas.marrSVessZv);
        // 3!

        // 4. формирование измерений палубных углов
        memcpy(BigMeas.marrMuWaveZv, trueMeasParams.marrMuWaveZv, 3 * sizeof(double));
        memcpy(BigMeas.marrMuZv, trueMeasParams.marrMuZv, 3 * sizeof(double));
        // 4!

        BigMeas.mTobr = mTObrabotki;
}
//--------------------------------------
// вычисление измеренного положения судна
// INPUT:
//arrSGskVessTrue[3]- вектор истинного положения судна в ГСК
//arrMu[3] - вектор палубных углов истинный
//arrMuZv[3] - вектор палубных углов измеренный
//OUTPUT:
// arrSGskVessZv[3] - выектор измеренного положения судна
void QImitMod:: imitateVesselPosition(double *arrSGskVessTrue,double *arrMu
                                      ,double *arrMuZv,double *arrSGskVessZv)
{
    // 1. вычиcление истинного положеиния GPS в ГСК
double arrTrueGPS_Position_GSK[3]={0.},arrTrueGPS_Position_KGSK[3]={0.}
                       , arrGPS_Position_GSK_Zv[3]={0.};
QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( mVess.mGps.marrPos, arrMu
                      , NULL,arrTrueGPS_Position_KGSK,3 );
MtrxSumMatrx(arrTrueGPS_Position_KGSK, arrSGskVessTrue,1, 3, arrTrueGPS_Position_GSK) ;

// ОТЛАДКА
//double matrPereh_PSK_V_KGSK[9] = {0.}, arr_kgsk[3] = {0.} , arr_gsk[3] = {0.};
//QCoordSystTrsf::calcMatr_PSK_v_KGSK_RightRot(arrMu, matrPereh_PSK_V_KGSK);
//trxMultMatrx( matrPereh_PSK_V_KGSK,3, 3, mVess.mGps.marrPos,1, arr_kgsk) ;
//MtrxSumMatrx(arr_kgsk, arrSGskVessTrue,1, 3, arr_gsk) ;

// !ОТЛАДКА
// 1!

// 2. формирование замера GPS
mVess.mGps.imitateMeasure(arrTrueGPS_Position_GSK, arrGPS_Position_GSK_Zv);
// 2!


  // 3. вычисление оценки положения GPS в КГСК на момент времени запросного сигнала
double arrPosGps_KGSK_Zv[3] = {0.};
QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( marrAprioriGpsPosParams, arrMuZv
                      , nullptr,arrPosGps_KGSK_Zv,3 );
// 3!

//4. формирование измерения положения судна на момент запросного сигнала
MtrxMinusMatrx(arrGPS_Position_GSK_Zv, arrPosGps_KGSK_Zv,1, 3, arrSGskVessZv);
// 4!
// ОТЛАДКА ЗАГЛУШКА
//mVess.mGps.imitateMeasure(arrSGskVessTrue, arrSGskVessZv);
// !ОТЛАДКА
}
//----------------------------------
void QImitMod:: createSeparatedPointsMeasures( const int QuantMeas, const double ZonaR_
                                               , QBigMeasure *arrMeas)
{
    for (int i =0;i < QuantMeas-1;++i)
    {
    // bear - пеленг корабля относительно маяка
    // double bear = -2. * M_PI/ (2.*((double)QuantMeas)) + 2. * M_PI/((double)QuantMeas) *((double)i);
    double bear =  2. * M_PI/((double)QuantMeas) *((double)i);
    // курсовой угол скорости корабля
    mVess.mQ0 = bear + M_PI/2.;

    mVess.marrVectSost[0] = marrAprioriBeaconPos[0] +  sin(bear) *ZonaR_ *0.4;
    mVess.marrVectSost[1] = marrAprioriBeaconPos[1] +  cos(bear) *ZonaR_ *0.4;
    mVess.marrVectSost[3]= 0.;
    mVess.marrVectSost[4]= 0.;
    mVess.marrVectSost[5]= 0.;
    mVess.mVQ = 0.;
    mVess.mVVess = 0.;
    mVess.mQ = mVess.mQ0; // угол курса
    mVess.mPsi = 0.; // угол килевой качки
    mVess.mTet = 0.; // угол бортовой качки
    mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];

    //memcpy(mVess.marrVectSost0, mVess.marrVectSost,9 * sizeof(double));
    mVess.tuneCurrentSinsMeasures();
    createSingleMeasure(arrMeas[i]);

    }

    // по центру

    mVess.mQ0 =  M_PI/2.;

    mVess.marrVectSost[0] = marrAprioriBeaconPos[0] ;
    mVess.marrVectSost[1] = marrAprioriBeaconPos[1] ;
    mVess.marrVectSost[3]= 0.;
    mVess.marrVectSost[4]= 0.;
    mVess.marrVectSost[5]= 0.;
    mVess.mVQ = 0.;
    mVess.mVVess = 0.;
    mVess.mQ = mVess.mQ0; // угол курса
    mVess.mPsi = 0.; // угол килевой качки
    mVess.mTet = 0.; // угол бортовой качки
    mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];

    //memcpy(mVess.marrVectSost0, mVess.marrVectSost,9 * sizeof(double));
    mVess.tuneCurrentSinsMeasures();
    createSingleMeasure(arrMeas[QuantMeas-1]);
}
//-----------------------------------------------------
void QImitMod::createQdrtTrajMeasures( const int QuantMeas, const double ZonaR_
                        , QBigMeasure *arrMeas)
{
    /*double *arr = new double [2 * QuantMeas];
    memset (arr, 0, 2 * QuantMeas * sizeof(double));
    int iarrParts[1] ={0};
    TURPolyLine pln(  1,QuantMeas,iarrParts,arr);*/
     const double VAlQdrtSideLength = ZonaR_* 0.4 *2.;//ZonaR_* 0.66 *2.;
    int numQdrtSide = 0;
    initPartOfQdrtTraj(numQdrtSide,marrAprioriBeaconPos, VAlQdrtSideLength);

    for (int i =0; i < QuantMeas; ++i)
    {
      int numQdrtSideCur = calcNumQdrtSide(marrAprioriBeaconPos,  VAlQdrtSideLength
                                           , numQdrtSide);
      if(numQdrtSideCur != numQdrtSide)
      {
        numQdrtSide = numQdrtSideCur ;
        initPartOfQdrtTraj(numQdrtSide,marrAprioriBeaconPos, VAlQdrtSideLength);
      }

      //pln.Points[i].X =mVess.marrVectSost[0];
      //pln.Points[i].Y =mVess.marrVectSost[1];
      createSingleMeasure(arrMeas[i]);
    }

}
//-----------------------------------------------
void QImitMod::createZigZagTrajMeasures( const int QuantMeas, const double ZonaR_
                        , QBigMeasure *arrMeas)
{
    /*double *arr = new double [2 * QuantMeas];
    memset (arr, 0, 2 * QuantMeas * sizeof(double));
    int iarrParts[1] ={0};
    TURPolyLine pln(  1,QuantMeas,iarrParts,arr);*/
     const double VAlQdrtSideLength = ZonaR_* 0.4 *2.;//ZonaR_* 0.66 *2.;
    int numZigZagPart = 0;
    initPartOfZigZagTraj(numZigZagPart,marrAprioriBeaconPos, VAlQdrtSideLength);

    for (int i =0; i < QuantMeas; ++i)
    {
      int numZigZagPartCur = calcNumPartOfZigZagTraj(marrAprioriBeaconPos,  VAlQdrtSideLength
                                           , numZigZagPart);
      if(numZigZagPartCur != numZigZagPart)
      {
        numZigZagPart = numZigZagPartCur ;
        initPartOfZigZagTraj(numZigZagPart,marrAprioriBeaconPos, VAlQdrtSideLength);
      }

      //pln.Points[i].X =mVess.marrVectSost[0];
      //pln.Points[i].Y =mVess.marrVectSost[1];
      createSingleMeasure(arrMeas[i]);
    }

}
//-------------------------------------------------
void  QImitMod::initPartOfZigZagTraj(const int numQdrtSide
                ,const double *arrAprioriBeaconPos,const double VAlQdrtSideLength)
{
    // bear - пеленг корабля относительно маяка
    double bear = - M_PI/ 4. +  M_PI/2. *((double)numQdrtSide);
    // курсовой угол скорости
    double course = 0.;
    switch(numQdrtSide % 8)
    {
    case 0:
         bear = - M_PI/ 4. ;
         course = M_PI/2.;
        break;
    case 1:
        bear =  M_PI/ 4. ;
        course = M_PI/ 4. + M_PI;
        break;
    case 2:
        bear =  M_PI/ 4. + M_PI;
        course = M_PI/2.;
        break;
    case 3:
        bear =  M_PI/ 4. + M_PI/2.;
        course = -M_PI/4.;
        break;
    case 4:
        bear = - M_PI/ 4. ;
        course = M_PI;
        break;
    case 5:
        bear =  M_PI/ 4. + M_PI;
        course = M_PI/4.;
        break;
    case 6:
        bear =  M_PI/ 4. ;
        course = M_PI;
        break;
    case 7:
        bear =  M_PI/ 4. + M_PI/2.;
        course = -M_PI/4.;

        break;

    default:
        break;
    }

// курсовой угол скорости корабля
mVess.mQ0 =course;
mVess.mQ = course; // угол курса

double hipoten = VAlQdrtSideLength/ 2.*sqrt(2.);

mVess.marrVectSost[0] = marrAprioriBeaconPos[0] +  sin(bear) *hipoten;
mVess.marrVectSost[1] = marrAprioriBeaconPos[1] +  cos(bear) *hipoten;
mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];
mVess.marrVectSost[3]= mVess.mVVess * sin(mVess.mQ0);
mVess.marrVectSost[4]= mVess.mVVess * cos(mVess.mQ0);
mVess.marrVectSost[5]= 0.;
mVess.mVQ = 0.;

mVess.tuneCurrentSinsMeasures();
}

//--------------------------------------------
int QImitMod::calcNumPartOfZigZagTraj(const double *arrQdrtCentre, const double VAlQdrtSideLength
                              , const int NUmQdrtSideCur)
{
    // координаты корабля в плоскости XY
    double valPntX = mVess.marrVectSost[0];
    double valPntY = mVess.marrVectSost[1];

    // координаты квадраты по оси OX
    double ax = arrQdrtCentre[0] - VAlQdrtSideLength/2;
    double bx = arrQdrtCentre[0] + VAlQdrtSideLength/2;

    // координаты квадраты по оси OY
    double ay = arrQdrtCentre[1] - VAlQdrtSideLength/2;
    double by = arrQdrtCentre[1] + VAlQdrtSideLength/2;

    if( (fabs(valPntX - arrQdrtCentre[0]) <= (VAlQdrtSideLength/2. + 0.1) )
     &&(fabs(valPntY - arrQdrtCentre[1]) <= (VAlQdrtSideLength/2. + 0.1) ))
    {
        return NUmQdrtSideCur;
    }

    return (NUmQdrtSideCur + 1) % 8 ;
}
 //-------------------------------------------------
void QImitMod::createDiamTrajMeasures( const int QuantMeas, const double ZonaR_
                        , QBigMeasure *arrMeas)
{
    /*double *arr = new double [2 * QuantMeas];
    memset (arr, 0, 2 * QuantMeas * sizeof(double));
    int iarrParts[1] ={0};
    TURPolyLine pln(  1,QuantMeas,iarrParts,arr);*/
     const double VAlDiamLength = ZonaR_* 0.4 *2.;//ZonaR_* 0.66 *2.;
    int numDiamPart = 0;
    initPartOfDiamTraj(numDiamPart,marrAprioriBeaconPos, VAlDiamLength);

    for (int i =0; i < QuantMeas; ++i)
    {
      double valHorizDist = sqrt( (mVess.marrVectSost[0] - marrAprioriBeaconPos[0]) * (mVess.marrVectSost[0] - marrAprioriBeaconPos[0])
              + (mVess.marrVectSost[1] - marrAprioriBeaconPos[1]) * (mVess.marrVectSost[1] - marrAprioriBeaconPos[1]));
       if (valHorizDist >= VAlDiamLength/2.)
       {
         numDiamPart++;
         initPartOfDiamTraj(numDiamPart,marrAprioriBeaconPos, VAlDiamLength);
       }

      createSingleMeasure(arrMeas[i]);
    }

}
//-------------------------------------------------
void  QImitMod::initPartOfDiamTraj(const int numDiamPart
                ,const double *arrAprioriBeaconPos,const double VAlDiamLength)
{
    // bear - пеленг корабля относительно маяка
    double bear =  M_PI/4. *((double)numDiamPart) + M_PI/8. * ((double)(numDiamPart%8));
    // курсовой угол скорости
    double course = M_PI + bear;

// курсовой угол скорости корабля
mVess.mQ0 =course;
mVess.mQ = course; // угол курса

mVess.marrVectSost[0] = marrAprioriBeaconPos[0] +  sin(bear) *VAlDiamLength / 2.;
mVess.marrVectSost[1] = marrAprioriBeaconPos[1] +  cos(bear) *VAlDiamLength / 2.;
mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];
mVess.marrVectSost[3]= mVess.mVVess * sin(mVess.mQ0);
mVess.marrVectSost[4]= mVess.mVVess * cos(mVess.mQ0);
mVess.marrVectSost[5]= 0.;
mVess.mVQ = 0.;

mVess.tuneCurrentSinsMeasures();
}
//-------------------------------------------------
void QImitMod::createLine6TrajMeasures( const int QuantMeas, const double ZonaR_
                       , QBigMeasure *arrMeas)
{
   /*double *arr = new double [2 * QuantMeas];
   memset (arr, 0, 2 * QuantMeas * sizeof(double));
   int iarrParts[1] ={0};
   TURPolyLine pln(  1,QuantMeas,iarrParts,arr);*/
    const double VAlDiamLength = ZonaR_* 0.4 ;//ZonaR_* 0.66 *2.;
   int numDiamPart = 0;
   initPartOfLine6Traj(numDiamPart,marrAprioriBeaconPos, VAlDiamLength);

   double val_part_X0 = mVess.marrVectSost[0];
   double val_part_Y0 = mVess.marrVectSost[1];

   for (int i =0; i < QuantMeas; ++i)
   {
       // ОТЛАДКА
       //if (39 == i)
       //{
       //    int uuu =0;
      // }
         // ОТЛАДКА !
     double valHorizDist = sqrt( (mVess.marrVectSost[0] - val_part_X0) * (mVess.marrVectSost[0] - val_part_X0)
             + (mVess.marrVectSost[1] - val_part_Y0) * (mVess.marrVectSost[1] - val_part_Y0));
      if (valHorizDist >= VAlDiamLength/ sqrt(2.))
      {
        numDiamPart++;
        initPartOfLine6Traj(numDiamPart,marrAprioriBeaconPos, VAlDiamLength);
        val_part_X0 = mVess.marrVectSost[0];
        val_part_Y0 = mVess.marrVectSost[1];
      }

     createSingleMeasure(arrMeas[i]);
   }
  int nnn =numDiamPart;
}
//-------------------------------------------------
void  QImitMod::initPartOfLine6Traj(const int numDiamPart
                ,const double *arrAprioriBeaconPos,const double VAlDiamLength)
{
    double arr_bear[] ={-M_PI/4., -M_PI/2.,-3./4.*M_PI
                         ,-M_PI/4.,  0., M_PI/4.

                        ,M_PI/4., M_PI/2.,3./4.*M_PI

                         , 3./4.*M_PI,M_PI, 5./ 4. *M_PI};

    double arr_Q[] ={M_PI/2., M_PI/2.,M_PI/2.
                         ,M_PI, M_PI, M_PI
                     ,-M_PI/2., -M_PI/2.,-M_PI/2.
                     ,0., 0.,0.
                       };

    double arr_coefdist []= {1., 1./sqrt(2.), 1.
                             ,1., 1./sqrt(2.), 1.

                            ,1., 1./sqrt(2.), 1.
                             ,1., 1./sqrt(2.), 1.};
    int inum = numDiamPart%12;
    // bear - пеленг корабля относительно маяка
    double bear =  arr_bear[inum];
    // курсовой угол скорости
    double course = arr_Q[inum];

    double coeff = arr_coefdist[inum];

// курсовой угол скорости корабля
mVess.mQ0 =course;
mVess.mQ = course; // угол курса

double valSdvig = VAlDiamLength/8.;
mVess.marrVectSost[0] = valSdvig + marrAprioriBeaconPos[0] +  sin(bear) *coeff *VAlDiamLength/2.;
mVess.marrVectSost[1] = valSdvig+ marrAprioriBeaconPos[1] +  cos(bear) *coeff *VAlDiamLength/2.;
mVess.marrVectSost[2] = -4. - mVess.mPlatform.marrPosParams[2];
mVess.marrVectSost[3]= mVess.mVVess * sin(mVess.mQ0);
mVess.marrVectSost[4]= mVess.mVVess * cos(mVess.mQ0);
mVess.marrVectSost[5]= 0.;
mVess.mVQ = 0.;

mVess.tuneCurrentSinsMeasures();
}

