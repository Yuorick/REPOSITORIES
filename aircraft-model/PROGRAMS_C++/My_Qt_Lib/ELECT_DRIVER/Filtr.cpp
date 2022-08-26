#include "Filtr.h"
#include <string.h>
#include "MeasurmentImitator.h"
#include "MatrixProccess.h"
#include <math.h>
//#include <float.h>


//extern double valOmReal;// отладка !!!

QFiltr::~QFiltr()
{
    if(marrCurEst) delete []marrCurEst ;
    marrCurEst = NULL ;
    if(marrCurK) delete []marrCurK ;
    marrCurK = NULL ;
    if(marrC) delete []marrC ;
    marrC = NULL ;

}
//---------------------------------------------------------------------------
QFiltr::QFiltr()
{
    mQntVars = 0;
    mQntObserv = 0 ;
    marrCurEst = NULL;
    marrCurK = NULL;
    marrC = NULL;
    mTCur = 0.;
    mh = 0.;
    mbInit = false;
}


//---------------------------------------------------------------------------

// конструктор копирования
 QFiltr::QFiltr (const  QFiltr &R)
 {
     mQntVars = R.mQntVars;
     mQntObserv = R.mQntObserv;
     mTCur = R.mTCur;
     mTCur = R.mTCur;
     mh = R.mh;
     mbInit = R.mbInit;

     if (R.marrCurEst != NULL)
     {
       if (!(marrCurEst = new double  [mQntVars]))
           abort();
     }

     if (R.marrCurK != NULL)
     {
       if (!(marrCurK = new double  [mQntVars * mQntVars]))
           abort();
     }

     if (R.marrC != NULL)
     {
       if (!(marrC = new double  [mQntVars * mQntObserv]))
           abort();
     }


     memcpy(marrCurEst, R.marrCurEst,  mQntVars * sizeof(double));
     memcpy(marrCurK, R.marrCurK,  mQntVars * mQntVars * sizeof(double));
     memcpy(marrC, R.marrC,  mQntVars * mQntObserv * sizeof(double));

 }

 // оператор присваивания
 QFiltr  &QFiltr::operator=( const QFiltr  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      mQntVars = R.mQntVars;
      mQntObserv = R.mQntObserv;
      mTCur = R.mTCur;
      mTCur = R.mTCur;
      mh = R.mh;
      mbInit = R.mbInit;

      if (R.marrCurEst != NULL)
      {
        if (!(marrCurEst = new double  [mQntVars]))
            abort();
      }

      if (R.marrCurK != NULL)
      {
        if (!(marrCurK = new double  [mQntVars * mQntVars]))
            abort();
      }

      if (R.marrC != NULL)
      {
        if (!(marrC = new double  [mQntVars * mQntObserv]))
            abort();
      }


      memcpy(marrCurEst, R.marrCurEst,  mQntVars * sizeof(double));
      memcpy(marrCurK, R.marrCurK,  mQntVars * mQntVars * sizeof(double));
      memcpy(marrC, R.marrC,  mQntVars * mQntObserv * sizeof(double));
     return *this ;
 }


 //-------------------------------------------------------------------

 // парам конструктор 1
 QFiltr:: QFiltr ( const int QntVars, const int QntObserv,const double *parrCurEst
  , const double *parrCurK, const double TCur, const double h
  , const bool bInit, const double *parrC)
 {
     mQntVars = QntVars;
     mQntObserv = QntObserv;
     mTCur = TCur;
     mh = h;
     mbInit = bInit;

     marrCurEst = NULL;
     marrCurK = NULL;
     marrC = NULL;


     if (parrCurEst != NULL)
     {
       if (!(marrCurEst = new double  [QntVars]))
           abort();
     }

     if (parrCurK != NULL)
     {
       if (!(marrCurK = new double  [QntVars * QntVars]))
           abort();
     }

     if (parrC != NULL)
     {
       if (!(marrC = new double  [QntVars * QntObserv]))
           abort();
     }

     memcpy(marrCurEst, parrCurEst,  QntVars * sizeof(double));
     memcpy(marrCurK, parrCurK,  QntVars * QntVars * sizeof(double));
     memcpy(marrC, parrC,  QntVars * QntObserv * sizeof(double));

 }
 //-------------------------------------------------------------------

 // парам конструктор 2
 QFiltr:: QFiltr ( const int QntVars, const int QntObserv,const double *parrCurEst
  , const double TCur, const double h  )
 {
     mQntVars = QntVars;
     mQntObserv = QntObserv;
     mTCur = TCur;
     mh = h;
     mbInit = true;;

     marrCurEst = NULL;
     marrCurK = NULL;
     marrC = NULL;


     if (parrCurEst != NULL)
     {
       if (!(marrCurEst = new double  [QntVars]))
           abort();
     }


     if (!(marrCurK = new double  [QntVars * QntVars]))
     {
         abort();
     }


       if (!(marrC = new double  [QntVars * QntObserv]))
       {
           abort();
       }


     memcpy(marrCurEst, parrCurEst,  QntVars * sizeof(double));
     memset(marrCurK, 0., QntVars *QntVars * sizeof(double));
     for (int i = 0; i < QntVars; ++i)
     {
       marrCurK[ i * QntVars + i ] = 0.00000000001;
     }

     memset(marrC, 0., QntVars * QntObserv * sizeof(double));
     for (int i = 0; i < QntObserv; ++i)
     {
       marrC[ i * QntVars + 1 + i] = 1.;
     }


 }
//--------------------------------
// инициализация
void QFiltr:: init_(const double valTcur,const double valh, const double*arrPhVect)
{
   mTCur =  valTcur;
   memcpy(marrCurEst, arrPhVect, mQntVars * sizeof(double));
   mh = valh;

   const double arrObsC[12] = {0., 1.,0.,0.
                        ,0.,0.,1.,0.
                        ,0.,0.,0.,1.};
   memcpy(marrC, arrObsC, 12 * sizeof(double));
   memset(marrCurK, 0, 16 * sizeof(double));
   for (int i = 0; i < 4; ++i)
   {
     marrCurK[ i * 4 + i ] = 0.00000000001;
   }
   mbInit = true;

}
/*
//------------------------------------------------
// обработка измерения ("учет замера", marrCurEst -экстраполированнок значение)
//arrA - матрица в правой части системы лин диф уравнений
void QFiltr::processMeasure(QMeasure MEasureCur,const double VAlTcur
                            ,double *arrA,double *arrBU,double *arrKww)
{
  // marrCurEst[0] = valOmReal ; // отладка !!!
  // memcpy(&(marrCurEst[1]), MEasureCur.marrYzv, 3 * sizeof(double));
    // матрица перехода
    double valH = VAlTcur - mTCur;
    double  *arr_L = new double[mQntVars * mQntVars];
   MatrxMultScalar(arrA, mQntVars, mQntVars, valH,arr_L);
   for(int i =0; i < mQntVars; ++i)
   {
     arr_L[ mQntVars * i + i] += 1.;
   }
   ///

   // экстраполяция фазового вектора
  // double arrPhVectExtr[4 ] = {0.},arrTemp0[4] = {0.}, arrBUt[4] = {0.};
  // MtrxMultMatrx(arr_L,4, 4,  marrCurEst,1, arrTemp0);


   //MatrxMultScalar(arrBU, 1, 4, valH,arrBUt);
  // MtrxSumMatrx(arrTemp0, arrBUt,4, 1, arrPhVectExtr) ;
   ///

   // экстраполяция коррел матрицы
   double *arrTemp1 = new double[mQntVars * mQntVars];
   double *arrTemp1_0 = new double[mQntVars * mQntVars];
   double *arrKExtr = new double[mQntVars * mQntVars];
   memcpy(arrTemp1_0, arr_L, mQntVars * mQntVars * sizeof(double));
   MtrxMultMatrx_MultMatrxTransp(arr_L,marrCurK,arrTemp1_0,mQntVars, arrTemp1);
   MtrxSumMatrx(arrTemp1, arrKww,mQntVars, mQntVars,  arrKExtr) ;
   ///

   // матрица коэффициентов усиления
   double* arrP = new double[mQntVars * mQntObserv];
   double *arrTemp2 = new double [mQntObserv * mQntObserv];
   double * arrTemp3 = new double [mQntObserv * mQntObserv];
   double * arrTemp3Inv = new double [mQntObserv * mQntObserv];
   MtrxMultMatrx_MultMatrxTransp(marrC,arrKExtr,marrC,mQntObserv,mQntVars, arrTemp2);
   MtrxSumMatrx(arrTemp2, MEasureCur.marrKYzv,mQntObserv, mQntObserv,  arrTemp3) ;

   for (int i = 0; i < 2; ++i)
   {
     arrTemp3[ i * mQntObserv + i] += 0.000001;
   }
   MatrxMultScalar(arrTemp3, mQntObserv, mQntObserv, 1.,arrTemp3);

   if (!InverseMtrx(arrTemp3, mQntObserv,arrTemp3Inv))
   {
       int iff =0;
   }

   MatrxMultScalar(arrTemp3Inv, mQntObserv, mQntObserv, 1.,arrTemp3Inv);
   double *arrT5 = new double [mQntVars * mQntObserv];
   MtrxMultMatrxTransp(arrKExtr,mQntVars, mQntVars, marrC,mQntObserv, arrT5) ;
   MtrxMultMatrx(arrT5,mQntVars, mQntObserv,  arrTemp3Inv,mQntObserv, arrP);
///

   // вычисление К(t|t)
   double *arrT6= new double[mQntVars * mQntVars];
   double *arrT7= new double[mQntVars * mQntVars];

   MtrxMultMatrx(arrP,mQntVars,mQntObserv,  marrC, mQntVars, arrT6);
   MtrxMultMatrx(arrT6,mQntVars, mQntVars,  arrKExtr,mQntVars, arrT7);
   MtrxMinusMatrx(arrKExtr, arrT7,mQntVars, mQntVars, marrCurK);

   ///

   // поправка оценки фазового вектора
    double *arrT9  = new double[mQntObserv];
    double *arrT10 = new double[mQntVars];
    double *arrDel = new double[mQntObserv];
    //

   MtrxMultMatrx(marrC,mQntObserv, mQntVars,  marrCurEst,1, arrT9);
   MtrxMinusMatrx(MEasureCur.marrYzv, arrT9,mQntObserv, 1, arrDel);
   MtrxMultMatrx(arrP,mQntVars, mQntObserv,  arrDel,1, arrT10);
   MtrxSumMatrx(marrCurEst, arrT10,1, mQntVars, marrCurEst) ;

   mTCur = VAlTcur;
   delete []arr_L;
   delete []arrTemp1;
   delete []arrTemp1_0;
   delete []arrKExtr;
   delete []arrP ;
   delete []arrTemp2;
   delete []arrTemp3 ;
   delete [] arrTemp3Inv;
   delete []arrT5;
   delete []arrT6;
   delete []arrT7;
   delete []arrT9;
   delete []arrT10;
   delete []arrDel;


}

*/
//------------------------------------------------
// обработка измерения ("учет замера", marrCurEst -экстраполированнок значение)
//arrA - матрица в правой части системы лин диф уравнений
void QFiltr::processMeasure(QDriverMeasure MEasureCur,const double VAlTcur
                            ,double *arrA,double *arrBU,double *arrKww)
{
  // marrCurEst[0] = valOmReal ; // отладка !!!
  // memcpy(&(marrCurEst[1]), MEasureCur.marrYzv, 3 * sizeof(double));
    // матрица перехода
    double valH = VAlTcur - mTCur;
    double  *arr_L = new double[mQntVars * mQntVars];
   MatrxMultScalar(arrA, mQntVars, mQntVars, valH,arr_L);
   for(int i =0; i < mQntVars; ++i)
   {
     arr_L[ mQntVars * i + i] += 1.;
   }
   ///

   // экстраполяция фазового вектора
  // double arrPhVectExtr[4 ] = {0.},arrTemp0[4] = {0.}, arrBUt[4] = {0.};
  // MtrxMultMatrx(arr_L,4, 4,  marrCurEst,1, arrTemp0);


   //MatrxMultScalar(arrBU, 1, 4, valH,arrBUt);
  // MtrxSumMatrx(arrTemp0, arrBUt,4, 1, arrPhVectExtr) ;
   ///

   // экстраполяция коррел матрицы
   double *arrTemp1 = new double[mQntVars * mQntVars];
   double *arrTemp1_0 = new double[mQntVars * mQntVars];
   double *arrKExtr = new double[mQntVars * mQntVars];
   memcpy(arrTemp1_0, arr_L, mQntVars * mQntVars * sizeof(double));
   MtrxMultMatrx_MultMatrxTransp(arr_L,marrCurK,arrTemp1_0,mQntVars, arrTemp1);
   MtrxSumMatrx(arrTemp1, arrKww,mQntVars, mQntVars,  arrKExtr) ;
   ///

   // матрица коэффициентов усиления
   double* arrP = new double[mQntVars * mQntObserv];
   double *arrTemp2 = new double [mQntObserv * mQntObserv];
   double * arrTemp3 = new double [mQntObserv * mQntObserv];
   double * arrTemp3Inv = new double [mQntObserv * mQntObserv];
   MtrxMultMatrx_MultMatrxTransp(marrC,arrKExtr,marrC,mQntObserv,mQntVars, arrTemp2);
   MtrxSumMatrx(arrTemp2, MEasureCur.marrKYzv,mQntObserv, mQntObserv,  arrTemp3) ;

   for (int i = 0; i < 2; ++i)
   {
     arrTemp3[ i * mQntObserv + i] += 0.000001;
   }
   MatrxMultScalar(arrTemp3, mQntObserv, mQntObserv, 1.,arrTemp3);

   if (!InverseMtrx(arrTemp3, mQntObserv,arrTemp3Inv))
   {
       int iff =0;
   }

   MatrxMultScalar(arrTemp3Inv, mQntObserv, mQntObserv, 1.,arrTemp3Inv);
   double *arrT5 = new double [mQntVars * mQntObserv];
   MtrxMultMatrxTransp(arrKExtr,mQntVars, mQntVars, marrC,mQntObserv, arrT5) ;
   MtrxMultMatrx(arrT5,mQntVars, mQntObserv,  arrTemp3Inv,mQntObserv, arrP);
///

   // вычисление К(t|t)
   double *arrT6= new double[mQntVars * mQntVars];
   double *arrT7= new double[mQntVars * mQntVars];

   MtrxMultMatrx(arrP,mQntVars,mQntObserv,  marrC, mQntVars, arrT6);
   MtrxMultMatrx(arrT6,mQntVars, mQntVars,  arrKExtr,mQntVars, arrT7);
   MtrxMinusMatrx(arrKExtr, arrT7,mQntVars, mQntVars, marrCurK);

   ///

   // поправка оценки фазового вектора
    double *arrT9  = new double[mQntObserv];
    double *arrT10 = new double[mQntVars];
    double *arrDel = new double[mQntObserv];
    //

   MtrxMultMatrx(marrC,mQntObserv, mQntVars,  marrCurEst,1, arrT9);
   MtrxMinusMatrx(MEasureCur.marrYzv, arrT9,mQntObserv, 1, arrDel);
   MtrxMultMatrx(arrP,mQntVars, mQntObserv,  arrDel,1, arrT10);
   MtrxSumMatrx(marrCurEst, arrT10,1, mQntVars, marrCurEst) ;

   mTCur = VAlTcur;
   delete []arr_L;
   delete []arrTemp1;
   delete []arrTemp1_0;
   delete []arrKExtr;
   delete []arrP ;
   delete []arrTemp2;
   delete []arrTemp3 ;
   delete [] arrTemp3Inv;
   delete []arrT5;
   delete []arrT6;
   delete []arrT7;
   delete []arrT9;
   delete []arrT10;
   delete []arrDel;


}

//------------------------------------------------
// обработка измерения
//arrL - матрица перехода
void QFiltr::processMeasure_Type2( double *ARrYZv,double *ARrK,const double VAlTcur
                            ,double *arrL,double *arrKww)
{


   // экстраполяция коррел матрицы
   double *arrTemp1 = new double[mQntVars * mQntVars];
   double *arrTemp1_0 = new double[mQntVars * mQntVars];
   double *arrKExtr = new double[mQntVars * mQntVars];
   memcpy(arrTemp1_0, arrL, mQntVars * mQntVars * sizeof(double));
   MtrxMultMatrx_MultMatrxTransp(arrL,marrCurK,arrTemp1_0,mQntVars, arrTemp1);
   MtrxSumMatrx(arrTemp1, arrKww,mQntVars, mQntVars,  arrKExtr) ;
   ///

   // матрица коэффициентов усиления
   double* arrP = new double[mQntVars * mQntObserv];
   double *arrTemp2 = new double [mQntObserv * mQntObserv];
   double * arrTemp3 = new double [mQntObserv * mQntObserv];
   double * arrTemp3Inv = new double [mQntObserv * mQntObserv];
   MtrxMultMatrx_MultMatrxTransp(marrC,arrKExtr,marrC,mQntObserv,mQntVars, arrTemp2);
   MtrxSumMatrx(arrTemp2, ARrK,mQntObserv, mQntObserv,  arrTemp3) ;



   for (int i = 0; i < mQntObserv; ++i)
   {
     arrTemp3[ i * mQntObserv + i] += 0.000001;
   }
   MatrxMultScalar(arrTemp3, mQntObserv, mQntObserv, 1.,arrTemp3);

   if (!InverseMtrx(arrTemp3, mQntObserv,arrTemp3Inv))
   {
       int iff =0;
   }

   MatrxMultScalar(arrTemp3Inv, mQntObserv, mQntObserv, 1.,arrTemp3Inv);
   double *arrT5 = new double [mQntVars * mQntObserv];
   MtrxMultMatrxTransp(arrKExtr,mQntVars, mQntVars, marrC,mQntObserv, arrT5) ;
   MtrxMultMatrx(arrT5,mQntVars, mQntObserv,  arrTemp3Inv,mQntObserv, arrP);
///

   // вычисление К(t|t)
   double *arrT6= new double[mQntVars * mQntVars];
   double *arrT7= new double[mQntVars * mQntVars];

   MtrxMultMatrx(arrP,mQntVars,mQntObserv,  marrC, mQntVars, arrT6);
   MtrxMultMatrx(arrT6,mQntVars, mQntVars,  arrKExtr,mQntVars, arrT7);
   MtrxMinusMatrx(arrKExtr, arrT7,mQntVars, mQntVars, marrCurK);

   ///

   // поправка оценки фазового вектора
    double *arrT9  = new double[mQntObserv];
    double *arrT10 = new double[mQntVars];
    double *arrDel = new double[mQntObserv];
    //

   MtrxMultMatrx(marrC,mQntObserv, mQntVars,  marrCurEst,1, arrT9);
   MtrxMinusMatrx(ARrYZv, arrT9,mQntObserv, 1, arrDel);
   MtrxMultMatrx(arrP,mQntVars, mQntObserv,  arrDel,1, arrT10);
   MtrxSumMatrx(marrCurEst, arrT10,1, mQntVars, marrCurEst) ;

   mTCur = VAlTcur;

   delete []arrTemp1;
   delete []arrTemp1_0;
   delete []arrKExtr;
   delete []arrP ;
   delete []arrTemp2;
   delete []arrTemp3 ;
   delete [] arrTemp3Inv;
   delete []arrT5;
   delete []arrT6;
   delete []arrT7;
   delete []arrT9;
   delete []arrT10;
   delete []arrDel;


}

