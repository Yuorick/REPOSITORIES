#include "PosSolver.h"

#include "BigMeasure.h"
#include "MatrixProccess.h"
#include <math.h>
#include "SubWaterBeam.h"
#include "LblSolver.h"
#include "Solver2D_Angs.h"
#include "Solver3D_Angs.h"
#include "CoordSystTrsf.h"

QPosSolver::~QPosSolver()
{

}
//---------------
QPosSolver::QPosSolver()
{
     mtblEstPrfl = TTable_1D();

    mDimX = 5;
    mDimY = 1;

    mDeepth=0.;
    mToler = 0.;

    mQuantIterMax = 0;

}
// конструктор копирования
 QPosSolver ::QPosSolver (const QPosSolver &R)
 {
     mVectBigMeasures.resize(R.mVectBigMeasures.size());
     mVectBigMeasures = R.mVectBigMeasures;
     mtblEstPrfl = R.mtblEstPrfl;     
     mDimX = R.mDimX;
     mDimY = R.mDimY;
     mDeepth = R.mDeepth;
     mToler = R.mToler;

     mQuantIterMax = R.mQuantIterMax;

 }


// оператор присваивания
QPosSolver &QPosSolver::operator=(const QPosSolver  &R)
{
    mVectBigMeasures.resize(R.mVectBigMeasures.size());
    mVectBigMeasures = R.mVectBigMeasures;
    mtblEstPrfl = R.mtblEstPrfl;   
    mDimX = R.mDimX;
    mDimY = R.mDimY;
    mDeepth = R.mDeepth;
    mToler = R.mToler;

    mQuantIterMax = R.mQuantIterMax;

  return *this ;
}

// парам констр 1
QPosSolver :: QPosSolver(const QBigMeasure *parrBigMeasures
        , const int  QntMeas    ,const TTable_1D tblEstPrfl
        , const double Deepth
        , const int DimX, const int DimY, const double Toler )

{
    mVectBigMeasures.resize(QntMeas);
    for (int i = 0; i < QntMeas; ++i)
    {
      mVectBigMeasures.replace(i, parrBigMeasures[i]);
    }

   mtblEstPrfl = tblEstPrfl;
   mDimX = DimX;
   mDimY = DimY;
   mDeepth = Deepth;
   mToler=Toler;
}


// парам констр 2
QPosSolver :: QPosSolver(QVector<QBigMeasure> VectBigMeasures
           ,const TTable_1D tblEstPrfl, const double Deepth
        , const int DimX, const int DimY, const double Toler )

{
    mVectBigMeasures.resize(VectBigMeasures.size());
    mVectBigMeasures = VectBigMeasures;


   mtblEstPrfl = tblEstPrfl;
   mDimX = DimX;
   mDimY = DimY;
   mDeepth = Deepth;
   mToler=Toler;
}


//--------------------------------------
// вычисление неизвестного вектора позиционирования методом наименьших квадратов
//arrX[0],arrX[1],arrX[2] - вектор параллакса антенны ПСК
//arrX[3],arrX[4],arrX[5] - вектор координат маяка в ГСК
//arrX[6],arrX[7],arrX[8] - вектор углов ориентации антенны в ПСК
// INPUT:
//arrX0[mDimX] - начальное приближение
// OUTPUT:
//arrXEst[mDimX] - решение
// valNeviaz0 - невязка
bool QPosSolver::estimateParams(  double *arrX0,double *arrXRez, double &valNeviaz0)
{
    memset(arrXRez, 0, 8 * sizeof(double));
    //0.
   double *arrMtrx_dFGr_po_dX_Inv = new double[mDimX * mDimX];
   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));


    // 1. формирование вектора начальгного приближения
    memcpy(arrXRez, arrX0, mDimX * sizeof(double));

    ///

    // 2. вычисление начальной невязки
    int quantGoodZamers0 = 0;
    valNeviaz0 = calcSumNeviazka(arrXRez, quantGoodZamers0);
    ///

    // 3. итреац процесс
   double *arrXCur = new double[mDimX];
   memcpy(arrXCur, arrXRez, mDimX * sizeof(double));

   double valNeviaz1=0.;
   double *arrDel = new double[mDimX];

   memset(arrDel, 0, mDimX* sizeof(double));

   bool breturn = false;
    for (int i =0; i < 400; ++i)
    {
     if(!doOneIteration(arrXRez, arrDel,arrMtrx_dFGr_po_dX_Inv))
     {
         delete []arrMtrx_dFGr_po_dX_Inv;
         delete []arrXCur;
         delete []arrDel;
         return false;
     }

     int len = mDimX;
     if (mDimX == 8)
     {
       len = 5;
     }
     double del1 = NormVect(arrDel, len);
     double del2 = (mDimX == 8)?NormVect(&arrDel[5], 3):0.;
        double coef = 0.2;//((i % 3)==0)?0.99:0.5;
        MatrxMultScalar( arrDel, mDimX, 1, coef, arrDel);
        MtrxMinusMatrx(arrXRez, arrDel,mDimX, 1, arrXCur);

        memcpy(arrXRez, arrXCur, mDimX * sizeof(double));
      if ((del1 < mToler)&&(del2 < 0.0005))
      {

          valNeviaz0 = valNeviaz1;
          breturn = true;
          break;
      }


        valNeviaz0 = valNeviaz1;
    }


    delete []arrMtrx_dFGr_po_dX_Inv;
    delete []arrXCur;
    delete []arrDel;
    return breturn;
}
//---------------------------------------------------------------
/*double QPosSolver::calcSumNeviazka(double *arrXCur, int &quantGoodZamers)
{
   double valSum = 0.;
   quantGoodZamers = 0;
   const int QntMeas = mVectBigMeasures.size();
   for (int i =0; i < QntMeas; ++i)
   {

       double valNeviazCur = 0.;
      if(!calcNeviazka(arrXCur,i, valNeviazCur))
      {
          continue;
      }
      quantGoodZamers++;
      valSum += valNeviazCur;
   }
    return valSum;
}
*/
double QPosSolver::calcSumNeviazka(double *arrXCur, int &quantGoodZamers)
{
   double valSum = 0.;
   quantGoodZamers = 0;
   const int QntMeas = mVectBigMeasures.size();
   double *parrNeviazCur = new double[QntMeas];
   bool *pbarrWellZamer = new bool[QntMeas];
   for (int i =0; i < QntMeas; ++i)
   {
      pbarrWellZamer[i] = true;
   }
   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
   #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
   { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


       for (int i =0; i < QntMeas; ++i)
       {


          if(!calcPartialNeviazkaSquare(arrXCur,i, parrNeviazCur[i]))
          {
              pbarrWellZamer[i] = false;
          }
          else
          {
           pbarrWellZamer[i] = true;
          }
       }

   } // OMP (начало блока, который выполняется в нескольких потоках !

   for (int i =0; i <QntMeas; ++i)
   {
     if(pbarrWellZamer[i])
     {
        valSum +=  parrNeviazCur[i];
        quantGoodZamers++;
     }
   }

   delete []parrNeviazCur;
   delete []pbarrWellZamer;
   return valSum;
}


double QPosSolver::calcArrNeviaz_Mean_and_Disp(double *arrXCur, int &quantGoodZamers
                                             , double *arrMean, double *arrDisp, double *arrNeviazSquare)
{
   double valSum = 0.;
   memset (arrDisp, 0, mDimY * sizeof(double));
   memset (arrNeviazSquare, 0, mDimY * sizeof(double));
   memset (arrMean, 0, mDimY * sizeof(double));
   quantGoodZamers = 0;
   const int QntMeas = mVectBigMeasures.size();
   double *parrNeviazCur = new double[QntMeas * mDimY];
   bool *pbarrWellZamer = new bool[QntMeas];
   for (int i =0; i < QntMeas; ++i)
   {
      pbarrWellZamer[i] = true;
   }
   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
   #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
   { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


       for (int i =0; i < QntMeas; ++i)
       {


          if(!calc_arrNeviaz(arrXCur, i, &parrNeviazCur[i * mDimY] ))
          {
              pbarrWellZamer[i] = false;
          }
          else
          {
           pbarrWellZamer[i] = true;
          }
       }

   } // OMP (начало блока, который выполняется в нескольких потоках !



   for (int i =0; i <QntMeas; ++i)
   {
     if(pbarrWellZamer[i])
     {
       MtrxSumMatrx(arrMean, &parrNeviazCur[i * mDimY],1, mDimY, arrMean) ;
       for (int j =0; j < mDimY; ++j)
       {
           arrNeviazSquare[j] += parrNeviazCur[i * mDimY +j] * parrNeviazCur[i * mDimY +j];
       }
        quantGoodZamers++;
     }
   }


   for (int i =0; i < mDimY; ++i)
   {
      arrMean[i] = arrMean[i]/ ((double)quantGoodZamers) ;
      arrNeviazSquare[i] = arrNeviazSquare[i]/ ((double)quantGoodZamers) ;
      arrDisp[i] = arrNeviazSquare[i] - arrMean[i] * arrMean[i];
   }


   for (int i =0; i < mDimY; ++i)
   {
      valSum += arrNeviazSquare[i];
   }

   delete []parrNeviazCur;
   delete []pbarrWellZamer;

   return valSum;
}



//---------------------------------------------------------
bool QPosSolver::calcPartialNeviazkaSquare(double *arrX,const int NUmBigMeasure, double &valNeviaz2)
{
  double *arrNeviaz = new double [mDimY];
  double *dArrNeviaz_po_dX = new double [mDimY * mDimX];
  if(!calc_arrNeviaz(arrX, NUmBigMeasure
                             , arrNeviaz ))
  {
      delete []arrNeviaz;
      return false;

  }


  valNeviaz2=  NormSquareVect(arrNeviaz, mDimY);
  return  true;
}
//-----------------------------------------------------------------
bool QPosSolver::doOneIteration(double *arrX0
           , double *arrDel, double *arrMtrx_dFGr_po_dX_Inv)
{

    double *arrFGr  = new double[mDimX];
    memset(arrFGr, 0,  mDimX * sizeof(double));
    double  *arrMtrx_dFGr_po_dX = new double[mDimX * mDimX];
    memset(arrMtrx_dFGr_po_dX, 0, mDimX * mDimX * sizeof(double));


    double valSumVeviaz2 = -1.;
  if (!calc_FGr_and_dFGr_po_dX(arrX0
                        , arrFGr, arrMtrx_dFGr_po_dX,valSumVeviaz2))
  {
      delete []arrFGr;
      delete []arrMtrx_dFGr_po_dX;
      return false;
  }

   double valScal = 1./3000.;
    MatrxMultScalar(arrMtrx_dFGr_po_dX, mDimX, mDimX,  valScal,arrMtrx_dFGr_po_dX);

    // PROBE
    double *arrPenalt = new double [mDimX];
    memset(arrPenalt, 0, mDimX * sizeof(double));
    double valPenalt = 100000.;
    calcArrPenalty(arrX0, valPenalt, arrPenalt);
    for (int  i =0; i < mDimX; ++i)
    {
        arrMtrx_dFGr_po_dX[i * mDimX +i] += arrPenalt[i];
    }
    delete []arrPenalt;
    // PROBE !
    bool brez =   InverseMtrx(arrMtrx_dFGr_po_dX, mDimX, arrMtrx_dFGr_po_dX_Inv);

    // ОТЛАДКА
  //  double arr[100]= {0.};
   // MtrxMultMatrx( arrMtrx_dFGr_po_dX,mDimX, mDimX, arrMtrx_dFGr_po_dX_Inv,mDimX, arr) ;

    // ! ОТЛАДКА   

    MatrxMultScalar(arrMtrx_dFGr_po_dX_Inv, mDimX, mDimX,  valScal,arrMtrx_dFGr_po_dX_Inv);

    if (!brez)
    {
        return false;
    }

    MtrxMultMatrx(arrMtrx_dFGr_po_dX_Inv,mDimX, mDimX, arrFGr,1, arrDel) ;

    delete []arrFGr;
    delete []arrMtrx_dFGr_po_dX;
    return true;

}

bool QPosSolver::calc_FGr_and_dFGr_po_dX(double *arrX0
               ,double * parrFGr,double * parrMtrx_dFGr_po_dX, double &valSumVeviaz2)
{
    memset(parrFGr, 0, mDimX * sizeof(double));
    memset(parrMtrx_dFGr_po_dX, 0, mDimX *mDimX * sizeof(double));


   valSumVeviaz2 = 0.;



   const int QntMeas = mVectBigMeasures.size();

   double *parrArr_fi = new double [QntMeas * mDimX];
   double *parrArr_df_po_dX = new double [QntMeas * mDimX* mDimX];
   double *parrArrNeviaz = new double [QntMeas * mDimY];
   bool *pbarrReturn = new bool [QntMeas];
   memset(pbarrReturn, true, QntMeas * sizeof(bool));
                   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
#pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


   for (int i = 0; i < QntMeas; ++i)
   {
       if (!calc_arrfi_and_arrdfi_po_dx(arrX0,i,&parrArrNeviaz[i * mDimY]
           , &parrArr_fi[i * mDimX],&parrArr_df_po_dX[i * mDimX *mDimX]))
       {

           pbarrReturn[i] =  false;
       }

   }


} // OMP (начало блока, который выполняется в нескольких потоках !

   for (int i =0; i < QntMeas; ++i)
   {
   if(!pbarrReturn[i])
   {
       delete []parrArr_fi;
       delete []parrArr_df_po_dX;
       delete []parrArrNeviaz;
       delete []pbarrReturn ;

       return false;
   }
   }

   for (int i = 0; i < QntMeas; ++i)
   {
       MtrxSumMatrx(parrFGr, &parrArr_fi[i * mDimX],1, mDimX, parrFGr) ;
       MtrxSumMatrx(parrMtrx_dFGr_po_dX, &parrArr_df_po_dX[i * mDimX *mDimX],mDimX, mDimX, parrMtrx_dFGr_po_dX) ;

       valSumVeviaz2 += NormSquareVect(&parrArrNeviaz[i * mDimY], mDimY);

   }

    delete []parrArr_fi;
    delete []parrArr_df_po_dX;
    delete []parrArrNeviaz;
   delete []pbarrReturn ;


    return true ;
}
//-------------------------------------------------
bool  QPosSolver::calc_arrfi_and_arrdfi_po_dx(double *arrX
                       ,const int NUmBigMeasure,double *arrNeviaz,double *arrfi,double *arrMtrx_dfi_po_dX)
{     

    // ОТЛАДКА
   // double arrAntPosParams[6] = {0.};
   // arrAntPosParams[0] = arrX[0];
   // arrAntPosParams[1] = arrX[1];
   // arrAntPosParams[2] = arrX[2];
    //arrAntPosParams[3] = arrX[5];
   // arrAntPosParams[4] = arrX[6];
    //arrAntPosParams[5] = arrX[7];
    //double arrGSK[3] = {0.};
   // double valTZv = -1.;
    //bool bb= transfMeasure_to_GSK(mtblEstPrfl,Meas
     //      , arrAntPosParams,NULL
     //      ,  arrGSK, &valTZv);

    // !ОТЛАДКА
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        return false;
    }
     // вектор приращений
    double arrPrirashen[] = {0.001,0.001,0.001,0.001,0.001,0.00001,0.00001,0.00001};

    double *arrXCur = new double [mDimX];
    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX * mDimX];
    double *arrNeviazCur = new double [mDimY];
    for (int i =0; i <mDimX; ++i)
    {
      memcpy(arrXCur, arrX, mDimX * sizeof(double));
      arrXCur[i]+= arrPrirashen[i];
      if (!calc_arrfi(arrXCur, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          delete []arrXCur;
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrTemp[mDimX * i]);
     // MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrMtrx_dfi_po_dX[mDimX * i]);
    }

    MatrTransp( arrTemp, mDimX, mDimX, arrMtrx_dfi_po_dX);

 delete []arrXCur;
 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;

    return true;

}
//-------------------------------------------------
bool  QPosSolver::calc_arrfi(double *arrX, const int NUmBigMeasure,double *arrNeviaz,double *arrfi)
{
    // 1. формирование вктора измерений и корреляц матрицы ошибок
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
   double *arrGrlsMeasure =new double [mDimY];
  double *arrD =new double [mDimY * mDimY];
   collectGrlsMeasure(Meas,arrGrlsMeasure,arrD);
    // !1

     // 2. вычисление  вектора невязок arrNeviaz[mDimY] и матрицы Якоби arr_dNeviaz_po_dX[mDimY * mDimX]
     // вектора невязок по X
     double *arr_dNeviaz_po_dX =new double [mDimY * mDimX];

     if (!calc_arrNeviaz_and_dArrNeviaz_po_dX(arrX, NUmBigMeasure
                                , arrNeviaz , arr_dNeviaz_po_dX))
     {
        delete []arrGrlsMeasure;
        delete []arrD;
         delete []arr_dNeviaz_po_dX;
         return false;
     }     
     // ОТЛАДКА
  //double arrNeviaz1[10] = {0.}, arrT[100]={0.},arr_dNeviaz_po_dX1_T[100] = {0.}, arr_dNeviaz_po_dX1[100] = {0.}
      //          ,arrX0_1[10] = {0.},arrdel[10] = {0.}
         //  , arrd[]= {0.001,0.001,0.001,0.001,0.001,0.00001,0.00001,0.00001};
   // for (int j =0; j < mDimX;++j)
  //  {
     //    memcpy(arrX0_1,arrX, mDimX * sizeof(double));
    //    arrX0_1[j] += arrd[j];
     //    calc_arrNeviaz_and_dArrNeviaz_po_dX(arrX0_1, NUmBigMeasure
        //                                 , arrNeviaz1 , arrT);

    // MtrxMinusMatrx(arrNeviaz1, arrNeviaz,1, mDimY, arrdel);
    //  MatrxMultScalar(arrdel, 1, mDimY, 1./arrd[j],&arr_dNeviaz_po_dX1_T[j * mDimY]);
   // }
   //  MatrTransp( arr_dNeviaz_po_dX1_T, mDimX, mDimY, arr_dNeviaz_po_dX1);
  // double ii = 0;
     // ОТЛАДКА !
      // 2!

     // 3. вычисление arrDInv
      double *arrDInv =new double [mDimY * mDimY];     
      InverseMtrx( arrD, mDimY, arrDInv) ;
     // fillE(arrD, mDimY);
      // 3!

       double *parrT = new double [mDimY];
       MtrxMultMatrx( arrNeviaz,1, mDimY, arrDInv,mDimY, parrT) ;
       MtrxMultMatrx( parrT,1, mDimY, arr_dNeviaz_po_dX,mDimX, arrfi) ;




  delete []arrGrlsMeasure;
  delete []arrD;
  delete []arr_dNeviaz_po_dX;
  delete []arrDInv;
  delete []parrT;
   return true;
}

//------------------------------------
void  QPosSolver::collectGrlsMeasure(const QBigMeasure &Meas
                          ,double *arrMeasure,double *arrD)
{

}
//-----------------------------------
bool QPosSolver::calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                                 , double *arrNeviaz , double *dArrNeviaz_po_dX)
{

}
//---------------------------------------------
/*
void  QPosSolver::calcNeviazArray(double * arrX0,double *parrBuffNeviaz,double *arrSKZ)
{
   double *parrf  = new double [mDimX] ;
   double *parr_df_po_dX = new double[mDimX * mDimX] ;
   const int QntMeas = mVectBigMeasures.size();
    for (int i = 0; i < QntMeas; ++i)
    {
        parrBuffNeviaz[4 * i] = ((double)i);
        if (!calc_arrfi_and_arrdfi_po_dx(arrX0, mVectBigMeasures.at(i),&parrBuffNeviaz[4 * i +1]
            , parrf,parr_df_po_dX))
        {
            parrBuffNeviaz[4 * i +1]=parrBuffNeviaz[4 * i +2]=parrBuffNeviaz[4 * i +3]= 1000.;
        }

    }

    delete []parrf;
    delete []parr_df_po_dX;

    double arrDisp[3] ={0.0000001,0.0000001,0.0000001};
    int countNBad = 0;
    for (int i =0; i < QntMeas; ++i)
    {
       for(int j =0; j < 3; ++j)
       {
           if(parrBuffNeviaz[4 * i + 1 + j] > 999.)
           {
              countNBad++;
              continue;
           }
           arrDisp[j] += parrBuffNeviaz[4 * i + 1 + j] * parrBuffNeviaz[4 * i + 1 + j];
       }
    }
    for(int j =0; j < 3; ++j)
    {
        arrSKZ[j] = sqrt(arrDisp[j]/((double)(QntMeas -countNBad -1.))) ;
    }
}
*/
//---------------------------------------------
void  QPosSolver::calcNeviazArray(double * arrX0,QVector <double> *vectBuffNeviaz,double *arrSKZ)
{
   double *parrf  = new double [mDimX] ;
   double *parr_df_po_dX = new double[mDimX * mDimX] ;
   const int QntMeas = mVectBigMeasures.size();
   (*vectBuffNeviaz).resize(4 *QntMeas);
   double arrNev[3] = {0.};
    for (int i = 0; i < QntMeas; ++i)
    {
       // parrBuffNeviaz[4 * i] = ((double)i);
        (*vectBuffNeviaz).replace(4 * i, ((double)i));

        if (!calc_arrfi_and_arrdfi_po_dx(arrX0, i,arrNev /*&parrBuffNeviaz[4 * i +1]*/
            , parrf,parr_df_po_dX))
        {
            //parrBuffNeviaz[4 * i +1]=parrBuffNeviaz[4 * i +2]=parrBuffNeviaz[4 * i +3]= 1000.;
            (*vectBuffNeviaz).replace(4 * i +1, 1000.);
            (*vectBuffNeviaz).replace(4 * i +2, 1000.);
            (*vectBuffNeviaz).replace(4 * i +3, 1000.);
        }
        else
        {
            (*vectBuffNeviaz).replace(4 * i +1, arrNev[0]);
            (*vectBuffNeviaz).replace(4 * i +2, arrNev[1]);
            (*vectBuffNeviaz).replace(4 * i +3, arrNev[2]);
        }

    }

    delete []parrf;
    delete []parr_df_po_dX;

    double arrDisp[3] ={0.0000001,0.0000001,0.0000001};
    int countNBad = 0;
    for (int i =0; i < QntMeas; ++i)
    {
        bool bZamerIsBad = false;
       for(int j =0; j < 3; ++j)
       {
          // if(parrBuffNeviaz[4 * i + 1 + j] > 999.)
          if(fabs((*vectBuffNeviaz).at(4 * i + 1 + j)) > 999.)
           {
              countNBad++;
              bZamerIsBad = true;
              break;
           }

           arrDisp[j] += (*vectBuffNeviaz).at(4 * i + 1 + j) * (*vectBuffNeviaz).at(4 * i + 1 + j);
       }
       if (bZamerIsBad)
       {
           continue;
       }
    }
    for(int j =0; j < 3; ++j)
    {
        arrSKZ[j] = sqrt(arrDisp[j]/((double)(QntMeas -countNBad -1.))) ;
    }
}
//-------------------------------------------
void QPosSolver::repaireZamerArray(double *arrX0)
{
    const int QntMeas = mVectBigMeasures.size();
   double arrDisp[3] ={0.000001,0.000001,0.000001}, arrSKZ[3] = {0.};
   //double *parrBuffNeviaz = new double [4 * QntMeas ];
   QVector <double>vectBuffNeviaz(4 * QntMeas );
   calcNeviazArray(arrX0,&vectBuffNeviaz,arrSKZ) ;


   //QVector<int>vctNumBad(mQntMeas);
   int quantBad = 0;
   bool bBad = false;
   for (  int ii =0; ii < QntMeas; ++ii)
   {
       bBad = false;
       int quantMeasCur =mVectBigMeasures.size();
        for (int i =0; i < quantMeasCur; ++i)
        {

           for(int j =0; j < 3; ++j)
           {
              //if (fabs(parrBuffNeviaz[4 * i + 1 + j]) > 2.5 *arrSKZ[j])
              if (fabs(vectBuffNeviaz.at(4 * i + 1 + j)) > 2.0 *arrSKZ[j])
              {
                  mVectBigMeasures.remove(i);
                  vectBuffNeviaz.remove(4 * i + 3);
                  vectBuffNeviaz.remove(4 * i + 2);
                  vectBuffNeviaz.remove(4 * i + 1);
                  vectBuffNeviaz.remove(4 * i);
                  bBad = true;
                  break;
              }
           }
           if(bBad)
           {
              break;
           }
        }

   }


   //delete []parrBuffNeviaz;

}
//--------------------------------------

void QPosSolver::unloadQvectorToDoubleArray(QVector <double> *vectInp,double *arrOut)
{
   for (int i =0; i < vectInp->size(); ++i)
   {
       arrOut[i] = vectInp->at(i);
   }
}
//------------------------------------------------

void QPosSolver::repaireZamerArray_(double *arrX0)
{
    const int QntMeas = mVectBigMeasures.size();
   double arrDisp[3] ={0.000001,0.000001,0.000001}, arrSKZ[3] = {0.};
   double *parrBuffNeviaz = new double [(1 + mDimY) * QntMeas ];
  // QVector <double>vectBuffNeviaz((1 + mDimY) * QntMeas );
   calcNeviazArray_(arrX0,parrBuffNeviaz,arrSKZ) ;


   //QVector<int>vctNumBad(mQntMeas);
   int quantBad = 0;
   bool bBad = false;
   for (  int ii = (QntMeas -1); ii >= 0; --ii)
   {



           for(int j =0; j < mDimY; ++j)
           {
              if (fabs(parrBuffNeviaz[(1 + mDimY) * ii + 1 + j]) > 2.5 *arrSKZ[j])
              //if (fabs(vectBuffNeviaz.at((1 + mDimY) + 1 + j)) > 2.5 *arrSKZ[j])
              {
                  mVectBigMeasures.remove(ii);
                  bBad = true;
                  break;
              }
           }

   }
delete []parrBuffNeviaz ;
}
//------------------------------------------------
//NUmVar  -номер переменной по которой идет чистка, нумерация от нуля
void QPosSolver::repaireZamerArray_(double *arrX0, const double VAlSigmaTreshold
                                    ,const int NUmVar, QVector<int> *pvctGoodZameresNum)
{

   const int QntMeas = mVectBigMeasures.size();
   double arrDisp[3] ={0.000001,0.000001,0.000001}, arrSKZ[3] = {0.};
   double *parrBuffNeviaz = new double [(1 + mDimY) * QntMeas ];

   calcNeviazArray_(arrX0,parrBuffNeviaz,arrSKZ) ;

   QVector <bool> bvctGoogZamer;
   bvctGoogZamer.resize(QntMeas);
   bvctGoogZamer.fill(true);
   int quantBad = 0;
   for (int i =0; i < QntMeas; ++i)
   {
       if (fabs(parrBuffNeviaz[(1 + mDimY)*i + 1 + NUmVar]) > VAlSigmaTreshold *arrSKZ[NUmVar])
       {
           bvctGoogZamer.replace(i,  false);
           quantBad ++;
           break;
       }

   }
   if (pvctGoodZameresNum != NULL)
   {
     pvctGoodZameresNum->resize(QntMeas - quantBad)  ;
   }
   int icur = 0;
   for (  int i =  (QntMeas-1); i >= 0; --i)
   {
     if (!bvctGoogZamer.at(i))
     {
        mVectBigMeasures.remove(i);
        continue;
     }

       if (pvctGoodZameresNum != NULL)
       {
         pvctGoodZameresNum->replace(icur, i) ;
         icur++;
       }

   }
  delete []parrBuffNeviaz;
}
//--------------------------------------
//parrBuffNeviaz[(1 + mDimY)*mVectBigMeasures.size()]
void  QPosSolver::calcNeviazArray_(double * arrX0, double *parrBuffNeviaz,double *arrSKZ)
{
   double *parrf  = new double [mDimX] ;
   double *parr_df_po_dX = new double[mDimX * mDimX] ;
   const int QntMeas = mVectBigMeasures.size();

   double *parrNev = new double [mDimY] ;
    for (int i = 0; i < QntMeas; ++i)
    {        
        parrBuffNeviaz[(1 + mDimY)* i] = (double)i;

        if (!calc_arrNeviaz(arrX0, i, parrNev ))
        {
            for (int k =0; k < mDimY; ++k)
            {
                parrBuffNeviaz[(1 + mDimY)* i +1 + k] = 1000.;
               //(*vectBuffNeviaz).replace((1 + mDimY)* i +1 + k, 1000.);
            }
        }
        else
        {
            for (int k =0; k < mDimY; ++k)
            {
                parrBuffNeviaz[(1 + mDimY)* i +1 + k] =parrNev[k];
              // (*vectBuffNeviaz).replace((1 + mDimY)* i +1 + k, parrNev[k]);
            }
        }

    }

    delete []parrf;
    delete []parr_df_po_dX;
    delete []parrNev;

    double arrDisp[3] ={0.0000001,0.0000001,0.0000001};
    int countNBad = 0;
    for (int i =0; i < QntMeas; ++i)
    {
        bool bZamerIsBad = false;
       for(int j =0; j < mDimY; ++j)
       {
           if(parrBuffNeviaz[(1 + mDimY)* i + 1 + j] > 999.)
         // if(fabs((*vectBuffNeviaz).at((1 + mDimY)* i + 1 + j)) > 999.)
           {
              countNBad++;
              bZamerIsBad = true;
              break;
           }

           arrDisp[j] += parrBuffNeviaz[(1 + mDimY)* i + 1 + j] * parrBuffNeviaz[(1 + mDimY)* i + 1 + j];
       }
       if (bZamerIsBad)
       {
           continue;
       }
    }
    for(int j =0; j < mDimY; ++j)
    {
        arrSKZ[j] = sqrt(arrDisp[j]/((double)(QntMeas -countNBad -1.))) ;
    }
}
//-------------------------
bool QPosSolver::calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                                                , double *arrNeviaz )
{
    return true;
}

//--------------------------------------------
void QPosSolver::fncRoughCleaningMeasures(double *arrX0,double valTreshold, const int NUmPeremennaya
                                          ,QVector<int>*pvctGoodZamersNum)
{
    double *parrNev = new double [mDimY];
    pvctGoodZamersNum->resize(mVectBigMeasures.size());
    int icur = 0;
    const int NUm0 = mVectBigMeasures.size();
    for (int i = (NUm0 -1); i >= 0; --i)
    {
        if (!calc_arrNeviaz(arrX0, i, parrNev ))
        {
            mVectBigMeasures.remove(i);
            continue;
        }

        if (fabs(parrNev[NUmPeremennaya]) > valTreshold)
        {
            mVectBigMeasures.remove(i);
            continue;
        }
        pvctGoodZamersNum->replace(icur,i);
        icur++;


    }
    pvctGoodZamersNum->resize(icur);
    delete []parrNev;
}
//---------------------------------
bool QPosSolver::calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(double *arrX0,double *arrFGr,double *arrMtrx_dFGr_po_dX_Inv)
{
    memset(arrFGr, 0,  mDimX * sizeof(double));
    double  *arrMtrx_dFGr_po_dX = new double[mDimX * mDimX];
    memset(arrMtrx_dFGr_po_dX, 0, mDimX * mDimX * sizeof(double));


    double valSumVeviaz2 = -1.;
  if (!calc_FGr_and_dFGr_po_dX(arrX0
                        , arrFGr, arrMtrx_dFGr_po_dX,valSumVeviaz2))
  {

      delete []arrMtrx_dFGr_po_dX;
      return false;
  }

    bool brez =   InverseMtrx(arrMtrx_dFGr_po_dX, mDimX, arrMtrx_dFGr_po_dX_Inv);

    if (!brez)
    {
        return false;
    }


    delete []arrMtrx_dFGr_po_dX;
    return true;
}
//-----------------------------------
void QPosSolver::createDblVect(const double *arr,const int len,QVector<double> &vct)
{
   vct.resize(len);
   for (int i=0; i< len; ++i)
   {
      vct.replace(i, arr[i]);
   }
}

//----------------------------------------
void QPosSolver:: createDblArrayFromDblVector(double *arr,QVector<double> vct)
{

   for (int i=0; i< vct.size(); ++i)
   {
      arr[i] = vct.at(i);
   }
}

//-------------------------------------
void QPosSolver::calcArrPenalty(const double *arrX0,const double valPenalt, double *arrPenalt)
{

}
//-------------------------------------
// задан порог по максимальной высоте антенны в ГСК (этоотрицательная величина!!!)
// требуется убрать все измерений у которых высота превосходит заданную
void QPosSolver::cleanZamersInAccordanceWithAntHeight(double *arrX0, const double VAlMaxHeightTreshold)
{
    const int QntMeas0 = mVectBigMeasures.size();
    double arrBeaconPos[3] = {0.};
    arrBeaconPos[0] = arrX0[3];
    arrBeaconPos[1] = arrX0[4];
    arrBeaconPos[2] = mDeepth;
    for(int j = QntMeas0 -1; j >= 0; --j)
    {

            QBigMeasure BigMeasure = mVectBigMeasures.at(j);
            // персчет вектора положения антенны в ГСК
            double arrAntY[3] = {0.};
            // создание матрицы перехода из   КГСК в ПСК
             double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
            MatrxMultScalar(BigMeasure.marrMuWaveZv, 1, 3, -1.,arrEilers);
            double matrPereh_PSK_V_KGSK[9] = {0} ;
            QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
            // вектор положения в ПСК-центр тяжести
            // вычисление вектора положениея в ПСК

            MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrX0, 1, arr_KGSK) ;

            MtrxSumMatrx(arr_KGSK, BigMeasure.marrSVessWaveZv,1, 3, arrAntY) ;
            // !
            if(arrAntY[2] >= VAlMaxHeightTreshold)
            {
               mVectBigMeasures.remove(j);
            }


            // создание матрицы перехода из   КГСК в ПСК

            MatrxMultScalar(BigMeasure.marrMuZv, 1, 3, -1.,arrEilers);

            QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
            // вектор положения в ПСК-центр тяжести
            // вычисление вектора положениея в ПСК

            MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrX0, 1, arr_KGSK) ;

            MtrxSumMatrx(arr_KGSK, BigMeasure.marrSVessZv,1, 3, arrAntY) ;
            // !
            if(arrAntY[2] >= VAlMaxHeightTreshold)
            {
               mVectBigMeasures.remove(j);
            }

    }
}
//-------------------------------------------------------------
bool QPosSolver::calc_darrfi_po_dSvess(double *arrX, const int NUmBigMeasure
                                       , const double del,double *parr_dfi_po_dSvess)
{
    double *arr_fi_lbl_0 = new double [ mDimX];

    double *arrNeviaz1 = new double[ mDimY];

    double *arr_fi_lbl_cur = new double [ mDimX];

    double *arr_dF_po_dP_T = new double [3 * mDimX];

    if(! calc_arrfi(arrX, NUmBigMeasure, arrNeviaz1,arr_fi_lbl_0))
    {
        delete []arr_fi_lbl_0;
        delete []arrNeviaz1;
        delete[]arr_fi_lbl_cur;
        delete[]arr_dF_po_dP_T;
        return false;
    }
    for (int j =0 ; j < 3; ++j)
    {
        QBigMeasure meas =  mVectBigMeasures.at(NUmBigMeasure);
        meas.marrSVessZv[j]  +=del;
         mVectBigMeasures.replace(NUmBigMeasure, meas);
         calc_arrfi(arrX,NUmBigMeasure, arrNeviaz1,arr_fi_lbl_cur);
        MtrxMinusMatrx(arr_fi_lbl_cur, arr_fi_lbl_0,1,  mDimX, arr_fi_lbl_cur);
        MatrxMultScalar(arr_fi_lbl_cur, 1,  mDimX, 1./del,&arr_dF_po_dP_T[j * mDimX]);

        meas.marrSVessZv[j]  -=del;
         mVectBigMeasures.replace(NUmBigMeasure, meas);
    }
    MatrTransp(arr_dF_po_dP_T,3,  mDimX, parr_dfi_po_dSvess);

    delete []arr_fi_lbl_0;
    delete []arrNeviaz1;
    delete[]arr_fi_lbl_cur;
    delete[]arr_dF_po_dP_T;
    return true;
}
//-------------------------------------

bool  QPosSolver::calc_arrfi_and_arrdfi_po_dt(double *arrX
                       ,const int NUmBigMeasure, const double valPrirashen,double *arrfi,double *arrdfi_po_dt)
{
   double arrNeviaz[8]= {0.};
    // !ОТЛАДКА
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        return false;
    }



    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX ];
    double arrNeviazCur[ 8] = {0.};


        QBigMeasure meas=mVectBigMeasures.at(NUmBigMeasure);
        meas.mTotvZv+= valPrirashen;
        mVectBigMeasures.replace(NUmBigMeasure, meas);

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {

          meas.mTotvZv -= valPrirashen;
          mVectBigMeasures.replace(NUmBigMeasure, meas);
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;

          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dt);
     meas.mTotvZv -= valPrirashen;
     mVectBigMeasures.replace(NUmBigMeasure, meas);
 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;


    return true;

}
//---------------------------------

bool  QPosSolver::calc_arrfi_and_arrdfi_po_dPsiTetta(double *arrX
                       ,const int NUmBigMeasure,double *arrPrirashen,double *arrdfi_po_dPsiTetta)
{
   double arrNeviaz[8]= {0.};
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrfi;
        return false;
    }




    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [2 * mDimX];
    double *arrNeviazCur = new double [mDimY];
    for (int i =0; i <2; ++i)
    {
        QBigMeasure meas= mVectBigMeasures.at(NUmBigMeasure);

        meas.marrMuZv[1 +i] += arrPrirashen[i];
        mVectBigMeasures.replace(NUmBigMeasure, meas);

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          meas.marrMuZv[1 +i] -= arrPrirashen[i];
          mVectBigMeasures.replace(NUmBigMeasure, meas);
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          delete []arrfi;

          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrTemp[mDimX * i]);
     meas.marrMuZv[1 +i] -= arrPrirashen[i];
     mVectBigMeasures.replace(NUmBigMeasure, meas);
    }

    MatrTransp( arrTemp, 2, mDimX, arrdfi_po_dPsiTetta);


 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrfi;


    return true;

}
//---------------------------------------------
bool  QPosSolver::calc_arrfi_and_arrdfi_po_dQ(double *arrX
                       ,const int NUmBigMeasure, const double valPrirashen,double *arrdfi_po_dQ)
{
   double arrNeviaz[8]= {0.};
   double *arrfi = new double[mDimX];
    // !ОТЛАДКА
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrfi;
        return false;
    }
    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX ];
    double *arrNeviazCur = new double [mDimY];


        QBigMeasure meas= mVectBigMeasures.at(NUmBigMeasure);
        meas.marrMuWaveZv[0] += valPrirashen;
        meas.marrMuZv[0] += valPrirashen;
        mVectBigMeasures.replace(NUmBigMeasure, meas);

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          meas.marrMuZv[0] -= valPrirashen;
          mVectBigMeasures.replace(NUmBigMeasure, meas);
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dQ);
     meas.marrMuZv[0] -= valPrirashen;
     mVectBigMeasures.replace(NUmBigMeasure, meas);
 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrfi;

    return true;

}
//-------------------------------------------------
bool  QPosSolver::calc_arrfi_and_arrdfi_po_dSgps(double *arrX
                       ,const int NUmBigMeasure,double *arrdfi_po_dSgps)
{
   double *arrNeviaz= new double [mDimY];
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrNeviaz;
        delete []arrfi;
        return false;
    }
     // вектор приращений
    double arrPrirashen[] = {0.1,0.1,0.1};


    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [3 * mDimX];
    double *arrNeviazCur = new double [mDimY];
    for (int i =0; i <3; ++i)
    {
        //QLblSolver LblSolverTemp = *this;
        QBigMeasure meas=mVectBigMeasures.at(NUmBigMeasure);

        meas.marrSVessZv[i]  += arrPrirashen[i];

        mVectBigMeasures.replace(NUmBigMeasure, meas);

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          meas.marrSVessZv[i]  -= arrPrirashen[i];
          mVectBigMeasures.replace(NUmBigMeasure, meas);
          delete []arrNeviaz;
          delete []arrfi;
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrTemp[mDimX * i]);
     meas.marrSVessZv[i]  -= arrPrirashen[i];
     mVectBigMeasures.replace(NUmBigMeasure, meas);
    }

    MatrTransp( arrTemp, 3, mDimX, arrdfi_po_dSgps);


 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrNeviaz;
 delete []arrfi;

    return true;

}

//-----------------------------------------

bool  QPosSolver::calc_arrfi_and_arrdfi_po_dAngSins(double *arrX
                       ,const int NUmBigMeasure,double *arrdfi_po_dAngSins)
{
   double *arrNeviaz = new double [mDimY];
   double *arrfi = new double [mDimX];
    // !ОТЛАДКА
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        return false;
    }
     // вектор приращений
    double arrPrirashen[] = {0.001,0.001,0.001};


    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX * mDimX];
    double *arrNeviazCur = new double [mDimY];
    for (int i =0; i <3; ++i)
    {
        QBigMeasure meas= mVectBigMeasures.at(NUmBigMeasure);
        meas.marrMuWaveZv[i] += arrPrirashen[i];
        meas.marrMuZv[i] += arrPrirashen[i];
        mVectBigMeasures.replace(NUmBigMeasure, meas);

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          meas.marrMuWaveZv[i] -= arrPrirashen[i];
          meas.marrMuZv[i] -= arrPrirashen[i];
          mVectBigMeasures.replace(NUmBigMeasure, meas);
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          delete []arrNeviaz;
          delete []arrfi;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrTemp[mDimX * i]);
     meas.marrMuWaveZv[i] -= arrPrirashen[i];
     meas.marrMuZv[i] -= arrPrirashen[i];
     mVectBigMeasures.replace(NUmBigMeasure, meas);

    }

    MatrTransp( arrTemp, 3, mDimX, arrdfi_po_dAngSins);


 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
    delete []arrNeviaz;
    delete []arrfi;


    return true;

}
//-----------------------------------------------------------------
bool  QPosSolver::calc_dFGr_po_dSins_(double* arrX0, double*  parrMtrx_dFGr_po_dSins)
{
    memset(parrMtrx_dFGr_po_dSins, 0, mDimX * 3 * sizeof(double));

    double * parrF0 = new double [mDimX];
    double * parrFcur = new double [mDimX];
    double * parrT = new double [mDimX];
    double * parrMtrx_dF_po_dX = new double [mDimX*mDimX];
    double valSumVeviaz2 = 0.;

    double*  parrMtrx_dFGr_po_dSins_T = new double [3*mDimX];

    calc_FGr_and_dFGr_po_dX(arrX0,parrF0,parrMtrx_dF_po_dX, valSumVeviaz2);

    double del = 0.0001;
    for (int i =0;i < 3; ++i)
    {
        QBigMeasure meas;
        for (int j = 0; j < mVectBigMeasures.size(); ++j)
        {
            meas = mVectBigMeasures.at(j);
            meas.marrMuZv[i] += del;
            meas.marrMuWaveZv[i] += del;
            mVectBigMeasures.replace(j,meas);
        }
        calc_FGr_and_dFGr_po_dX(arrX0,parrFcur, parrMtrx_dF_po_dX, valSumVeviaz2);
        MtrxMinusMatrx(parrFcur, parrF0,1, mDimX, parrT);
        MatrxMultScalar(parrT, 1, mDimX, 1./del,&parrMtrx_dFGr_po_dSins_T[i * mDimX]);

        for (int j = 0; j < mVectBigMeasures.size(); ++j)
        {
        meas = mVectBigMeasures.at(j);
        meas.marrMuZv[i] -= del;
        meas.marrMuWaveZv[i] -= del;
        mVectBigMeasures.replace(j,meas);
        }
    }

    MatrTransp( parrMtrx_dFGr_po_dSins_T, 3, mDimX, parrMtrx_dFGr_po_dSins);

    delete []parrF0 ;
    delete []parrFcur;
    delete []parrMtrx_dF_po_dX;
    delete []parrMtrx_dFGr_po_dSins_T;
    delete []parrT;
    return true ;

}
//-----------------------------------------------------------------
bool  QPosSolver::calc_dFGr_po_dProfile(double* arrX0, double* parrMtrx_dFGr_po_dProfile)
{
    memset(parrMtrx_dFGr_po_dProfile, 0, mDimX  * sizeof(double));

    double * parrF0 = new double [mDimX];
    double * parrFcur = new double [mDimX];
    double * parrT = new double [mDimX];
    double * parrMtrx_dF_po_dX = new double [mDimX*mDimX];
    double valSumVeviaz2 = 0.;

    calc_FGr_and_dFGr_po_dX(arrX0,parrF0,parrMtrx_dF_po_dX, valSumVeviaz2);

    double del = 0.1;


        for (int j = 0; j < mtblEstPrfl.mNumCols; ++j)
        {
            mtblEstPrfl.mparrVal[j] += del;
        }
        calc_FGr_and_dFGr_po_dX(arrX0,parrFcur, parrMtrx_dF_po_dX, valSumVeviaz2);
        MtrxMinusMatrx(parrFcur, parrF0,1, mDimX, parrT);
        MatrxMultScalar(parrT, 1, mDimX, 1./del,parrMtrx_dFGr_po_dProfile);

        for (int j = 0; j < mtblEstPrfl.mNumCols; ++j)
        {
            mtblEstPrfl.mparrVal[j] -= del;
        }



    delete []parrF0 ;
    delete []parrFcur;
    delete []parrMtrx_dF_po_dX;   
    delete []parrT;
    return true ;
}
