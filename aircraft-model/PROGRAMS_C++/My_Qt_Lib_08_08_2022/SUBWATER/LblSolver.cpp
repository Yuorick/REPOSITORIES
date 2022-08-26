#include "LblSolver.h"
#include "BigMeasure.h"
#include "SubWaterBeam.h"
#include "MatrixProccess.h"
//#include "DataExchange.h"

QLblSolver::QLblSolver():QPosSolver()
{
    mDimX = 5;
    mDimY =1;
}

// конструктор копирования
QLblSolver :: QLblSolver (const  QLblSolver &R):QPosSolver( R)
 {


 }

 // оператор присваивания
QLblSolver  &QLblSolver::operator=( const QLblSolver  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QPosSolver:: operator= (R);


     return *this ;
 }

// парам конструктор 1
QLblSolver:: QLblSolver (const QBigMeasure *parrBigMeasures
                         , const int  QntMeas    ,const TTable_1D tblEstPrfl
                         , const double Deepth
                         , const double Toler)
:QPosSolver (parrBigMeasures,  QntMeas    ,tblEstPrfl
             , Deepth, 5,  1,Toler)
{

}
//------------------------------------------

QLblSolver:: QLblSolver(QVector<QBigMeasure> VectBigMeasures
           ,const TTable_1D tblEstPrfl, const double Deepth
        , const double Toler )
    :QPosSolver(VectBigMeasures,tblEstPrfl,  Deepth
             , 5,  1, Toler )
{

}


void  QLblSolver::collectGrlsMeasure(const QBigMeasure &Meas
                                             ,double *arrMeasure,double *arrD)
{
   arrMeasure[0] = Meas.mTotvZv -  Meas.mTzaprZv - Meas.mTobr;
   arrD[0] = Meas.mSig_t*Meas.mSig_t;
}
//----------------------------------
//задан вектор неизвестных параметров arrX
//структура текущего измерения
// требуется вычислить какое будет измерение антенны
//INPUT:
//arrX[mDimX] -параметры
//Meas - замер
// OUTPUT:
//arrf[mDimX]-измерение
bool QLblSolver::calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                           , double *arrNeviaz , double *dArrNeviaz_po_dX)

{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
    // 1. вычисление времени запросного сигнала и его градиента
    double val_tZapr = -1.,  arrdtZapr_po_dX[6] = {0.};
    // выектор положения маяка
    double arrPosBeacon[3]={0.};
    arrPosBeacon[0] = arrX[3];
    arrPosBeacon[1] = arrX[4];
    arrPosBeacon[2] = mDeepth;

    if(! calc_t_and_dt_po_dX(mtblEstPrfl,  Meas.marrMuWaveZv, arrX
                             ,Meas.marrSVessWaveZv,arrPosBeacon
                                      ,val_tZapr, arrdtZapr_po_dX))
    {
        return false;
    }

    // 2. вычисление времени ответного сигнала и его градиента
    double val_tOtv = -1.,  arrdtOtv_po_dX[6] = {0.};

    if(! calc_t_and_dt_po_dX(mtblEstPrfl,  Meas.marrMuZv, arrX
                             ,Meas.marrSVessZv,arrPosBeacon
                                      ,val_tOtv, arrdtOtv_po_dX))
    {
        return false;
    }


    MtrxSumMatrx(arrdtZapr_po_dX, arrdtOtv_po_dX,1, 5, dArrNeviaz_po_dX) ;

    double coeff = 1.;
    MatrxMultScalar(dArrNeviaz_po_dX, 1, 5, coeff,dArrNeviaz_po_dX);
    arrNeviaz[0] =  coeff *(Meas.mTobr + val_tOtv + val_tZapr - (Meas.mTotvZv - Meas.mTzaprZv ));




    calc_2D_arrNeviaz_and_dArrNeviaz_po_dX(arrX, Meas
                               , arrNeviaz , dArrNeviaz_po_dX);

    return true;

}
bool QLblSolver::calc_2D_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, QBigMeasure Meas,
                            double *arrNeviaz , double *dArrNeviaz_po_dX)
{
    return true;
}
//------------------------------------
//----------------------------------
//задан вектор неизвестных параметров arrX
//структура текущего измерения
// требуется вычислить какое будет измерение антенны
//INPUT:
//arrX[mDimX] -параметры
//Meas - замер
// OUTPUT:
//arrf[mDimX]-измерение
bool QLblSolver::calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                           , double *arrNeviaz )

{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
    // 1. вычисление времени запросного сигнала и его градиента
    double val_tZapr = -1.,  arrdtZapr_po_dX[6] = {0.};
    // выектор положения маяка
    double arrPosBeacon[3]={0.};
    arrPosBeacon[0] = arrX[3];
    arrPosBeacon[1] = arrX[4];
    arrPosBeacon[2] = mDeepth;

    if(!calc_t_for_Vess_gsk(mtblEstPrfl,  Meas.marrMuWaveZv,  arrX
                            , Meas.marrSVessWaveZv,arrPosBeacon
                                     ,val_tZapr) )
    {
        return false;
    }

    // 2. вычисление времени ответного сигнала и его градиента
    double val_tOtv = -1.;
    if(!calc_t_for_Vess_gsk(mtblEstPrfl,  Meas.marrMuZv,  arrX
                            , Meas.marrSVessZv,arrPosBeacon
                                     ,val_tOtv) )
    {
        return false;
    }

    arrNeviaz[0] =  (Meas.mTobr + val_tOtv + val_tZapr - (Meas.mTotvZv - Meas.mTzaprZv ));

    return true;

}
//---------------------------------
void QLblSolver::calc_dX_po_dSgps(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dSgps)
{
    double  *parrMtrx_dFGr_po_dSgps = new double[mDimX *3];
    memset (parrMtrx_dFGr_po_dSgps, 0, mDimX *3 * sizeof(double));
    calc_FGr_and_dFGr_po_dSgps(arrX0, parrMtrx_dFGr_po_dSgps);
    MtrxMultMatrx( arrMtrx_dFGr_po_dX_Inv,mDimX, mDimX, parrMtrx_dFGr_po_dSgps,3, arr_dX_po_dSgps) ;
    MatrxMultScalar(arr_dX_po_dSgps, 3, mDimX, -1.,arr_dX_po_dSgps);
    delete []parrMtrx_dFGr_po_dSgps;

}
//-------------------------------------------
bool QLblSolver::calc_FGr_and_dFGr_po_dSgps(double *arrX0,double * parrMtrx_dFGr_po_dSgps)
{
   memset(parrMtrx_dFGr_po_dSgps, 0, mDimX *3 * sizeof(double));

   const int QntMeas = mVectBigMeasures.size();


   double *parrArr_df_po_dSgps = new double [QntMeas * mDimX* 3];

   bool *barRreturn = new bool[QntMeas];
   memset(barRreturn, false,QntMeas * sizeof(bool));
                   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
  #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


   for (int i = 0; i < QntMeas; ++i)
   {


       if (!calc_arrfi_and_arrdfi_po_dSgps(arrX0
                                           ,i,&parrArr_df_po_dSgps[i * mDimX* 3])
              )
       {

           barRreturn[i] =  true;
       }


   }


} // OMP (начало блока, который выполняется в нескольких потоках !
   for (int i =0; i < QntMeas; ++i)
   {
       if(barRreturn[i])
       {
           delete []parrArr_df_po_dSgps;
           delete []barRreturn;
           return false;
       }
   }
   for (int i = 0; i < QntMeas; ++i)
   {

       MtrxSumMatrx(parrMtrx_dFGr_po_dSgps, &parrArr_df_po_dSgps[i * mDimX* 3],mDimX, 3, parrMtrx_dFGr_po_dSgps) ;



   }


    delete []parrArr_df_po_dSgps;
   delete []barRreturn;

    return true ;
}


//-----------------------------------------------
void QLblSolver::calc_dX_po_dHBeacon(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dHBeacon)
{
    double  parrMtrx_dFGr_po_dHBeacon[5]= {0.};
    calc_dFGr_po_dHBeacon(arrX0, parrMtrx_dFGr_po_dHBeacon);
    MtrxMultMatrx( arrMtrx_dFGr_po_dX_Inv,mDimX, mDimX, parrMtrx_dFGr_po_dHBeacon,1, arr_dX_po_dHBeacon) ;
    MatrxMultScalar(arr_dX_po_dHBeacon, 1, mDimX, -1.,arr_dX_po_dHBeacon);
}

//-------------------------------------------
bool QLblSolver::calc_dFGr_po_dHBeacon(double *arrX0,double * parrMtrx_dFGr_po_dHBeacon)
{

    memset(parrMtrx_dFGr_po_dHBeacon, 0, mDimX  * sizeof(double));

   const int QntMeas = mVectBigMeasures.size();


   double *parrArr_df_po_dHBeacon = new double [QntMeas * mDimX];

   bool *barRreturn = new bool[QntMeas];
   memset(barRreturn, false,QntMeas * sizeof(bool));
                   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
  #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


   for (int i = 0; i < QntMeas; ++i)
   {


       if (!calc_arrfi_and_arrdfi_po_dHBeacon(arrX0
                                           ,i,&parrArr_df_po_dHBeacon[i * mDimX])
              )
       {

           barRreturn[i] = true;
       }

   }


} // OMP (начало блока, который выполняется в нескольких потоках !
   for (int i =0; i < QntMeas; ++i)
   {
       if(barRreturn[i])
       {
           delete []parrArr_df_po_dHBeacon;
           delete []barRreturn;
           return false;
       }
   }
   for (int i = 0; i < QntMeas; ++i)
   {

       MtrxSumMatrx(parrMtrx_dFGr_po_dHBeacon, &parrArr_df_po_dHBeacon[i * mDimX],mDimX, 1, parrMtrx_dFGr_po_dHBeacon) ;



   }


   delete []parrArr_df_po_dHBeacon;
   delete []barRreturn;

    return true ;
}
//-------------------------------------------------
bool  QLblSolver::calc_arrfi_and_arrdfi_po_dHBeacon(double *arrX
                       ,const int NUmBigMeasure,double *arrdfi_po_dHBeacon)
{
   double *arrNeviaz = new double [mDimY];
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrNeviaz;
        delete []arrfi ;
        return false;
    }
     // вектор приращений
    double valPrirashen = 0.01;


    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];    
    double *arrNeviazCur = new double [mDimY];


        QLblSolver LblSolverTemp = *this;

        LblSolverTemp.mDeepth += valPrirashen;


      if (!LblSolverTemp.calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          delete []arrNeviaz;
          delete []arrfi ;
          delete []arrfiCur;
          delete []arrRazn;         
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dHBeacon);

 delete []arrfiCur;
 delete []arrRazn;
 delete []arrNeviazCur;
 delete []arrNeviaz;
 delete []arrfi ;

    return true;

}
//-------------------------------------------------
// Вычислениекореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
//
//
bool  QLblSolver::calc_Kgps(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kgps_per_1_m)
{

    memset(arr_Kgps_per_1_m, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();


    double *parrArr_df_po_dSgps = new double [QntMeas * mDimX*3];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrfi_and_arrdfi_po_dSgps(arrX0
                                            ,i,&parrArr_df_po_dSgps[i * mDimX* 3])
               )
        {

            barRreturn[i] =  true;
        }

    }


 } // OMP (начало блока, который выполняется в нескольких потоках !
    for (int i =0; i < QntMeas; ++i)
    {
        if(barRreturn[i])
        {
            delete []parrArr_df_po_dSgps;
            delete []barRreturn;
            return false;
        }
    }

    double *arrTemp0 = new double [mDimX *mDimX];
    double *arrTemp1 = new double [mDimX *mDimX];
    memset(arrTemp1, 0, sizeof(double) *mDimX *mDimX);

    for (int i = 0; i < QntMeas; ++i)
    {

        calcF_Mult_FTransp(&parrArr_df_po_dSgps[i * mDimX* 3],mDimX, 3, arrTemp0);
        MtrxSumMatrx(arrTemp1, arrTemp0,mDimX, mDimX, arrTemp1) ;

    }


    delete []parrArr_df_po_dSgps;
    delete []barRreturn;


    calcF_D_FTransp(arrMtrx_dFGr_po_dX_Inv,arrTemp1,mDimX, arr_Kgps_per_1_m);
    delete []arrTemp0;
    delete []arrTemp1;

     return true ;
}
//--------------------------------

// Вычислениекореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
bool  QLblSolver::calc_Kt(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kt)
{
    memset(arr_Kt, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();


    double *parrArr_df_po_dt = new double [QntMeas * mDimX];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {

        double *parrArr_fi = new double [ mDimX];
        if (!calc_arrfi_and_arrdfi_po_dt(arrX0
                                            ,i,0.1, parrArr_fi,&parrArr_df_po_dt[i * mDimX])
               )
        {

             barRreturn[i] =  true;
        }
       delete []parrArr_fi;
    }


 } // OMP (начало блока, который выполняется в нескольких потоках !
    for (int i =0; i < QntMeas; ++i)
    {
        if(barRreturn[i])
        {

            delete []parrArr_df_po_dt;
            delete []barRreturn;
            return false;
        }
    }


    double *arrTemp0 = new double [mDimX *mDimX];
    double *arrTemp1 = new double [mDimX *mDimX];
    memset(arrTemp1, 0, sizeof(double) *mDimX *mDimX);
    double *arrTemp2 = new double [mDimX *mDimX];

    for (int i = 0; i < QntMeas; ++i)
    {
        calcF_Mult_FTransp(&parrArr_df_po_dt[i * mDimX],mDimX, 1, arrTemp0);
        MatrxMultScalar(arrTemp0, mDimX, mDimX, mVectBigMeasures.at(i).mSig_t *  mVectBigMeasures.at(i).mSig_t,arrTemp0);
        MtrxSumMatrx(arrTemp1, arrTemp0,mDimX, mDimX, arrTemp1) ;

    }


    delete []parrArr_df_po_dt;
    delete []barRreturn;


    calcF_D_FTransp(arrMtrx_dFGr_po_dX_Inv,arrTemp1,mDimX, arr_Kt);
    delete []arrTemp0;
    delete []arrTemp1;
    delete []arrTemp2;

     return true ;
}

//---------------------------------------

bool  QLblSolver::calc_dX_po_dAngSins(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dAngSins)
{
    memset(arr_dX_po_dAngSins, 0, mDimX *3 * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();


    double *parrArr_df_po_dAngSins = new double [QntMeas * mDimX* 3];
    double * parrMtrx_dFGr_po_dAngSins = new double [mDimX* 3];
    memset (parrMtrx_dFGr_po_dAngSins, 0, mDimX* 3 * sizeof(double));

    bool *barRreturn = new bool [QntMeas];
    memset (barRreturn , false, QntMeas *sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrfi_and_arrdfi_po_dAngSins(arrX0
             ,i,&parrArr_df_po_dAngSins[i * mDimX* 3])
               )
        {

            barRreturn [i] =   true;
        }

    }


 } // OMP (начало блока, который выполняется в нескольких потоках !

    for (int i = 0; i < QntMeas; ++i)
    {
         if (barRreturn [i])
         {

        delete []parrArr_df_po_dAngSins;
        delete []parrMtrx_dFGr_po_dAngSins;
        delete []barRreturn;

        return false;
        }
    }

    for (int i = 0; i < QntMeas; ++i)
    {

        MtrxSumMatrx(parrMtrx_dFGr_po_dAngSins, &parrArr_df_po_dAngSins[i * mDimX* 3],mDimX, 3, parrMtrx_dFGr_po_dAngSins) ;
    }


     delete []parrArr_df_po_dAngSins;
     delete []barRreturn;

    // ОТЛАДКА
   // double parrMtrx_dFGr_po_dSins[15] = {0.};
   // calc_dFGr_po_dSins_(arrX0,parrMtrx_dFGr_po_dSins);
    // !ОТЛАДКА
    MtrxMultMatrx( arrMtrx_dFGr_po_dX_Inv,mDimX, mDimX, parrMtrx_dFGr_po_dAngSins,3, arr_dX_po_dAngSins) ;
    MatrxMultScalar(arr_dX_po_dAngSins, 3, mDimX, -1.,arr_dX_po_dAngSins);

    delete []parrMtrx_dFGr_po_dAngSins;

     return true ;
}

//-------------------------------------------------
// Вычисление кореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
//
//
bool  QLblSolver::calc_KsinsQ(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KsinsQ_per_1_rad)
{

    memset(arr_KsinsQ_per_1_rad, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();

   // double *parrArr_fi = new double [QntMeas * mDimX];
    double *parrArr_df_po_dTetta = new double [QntMeas * mDimX];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrfi_and_arrdfi_po_dQ(arrX0
                                            ,i, 0.001,&parrArr_df_po_dTetta[i * mDimX])
               )
        {

            barRreturn[i] =  true;
        }

    }


} // OMP (начало блока, который выполняется в нескольких потоках !
    for (int i = 0; i < QntMeas; ++i)
    {
         if (barRreturn [i])
         {


             delete []parrArr_df_po_dTetta;
             delete []barRreturn;

        return false;
        }
    }

    double *arrTemp0 = new double [mDimX *mDimX];
    double *arrTemp1 = new double [mDimX *mDimX];
    memset(arrTemp1, 0, sizeof(double) *mDimX *mDimX);
    double *arrTemp2 = new double [mDimX *mDimX];

    for (int i = 0; i < QntMeas; ++i)
    {
        calcF_Mult_FTransp(&parrArr_df_po_dTetta[i * mDimX],mDimX, 1, arrTemp0);
        MtrxSumMatrx(arrTemp1, arrTemp0,mDimX, mDimX, arrTemp1) ;

    }


    delete []parrArr_df_po_dTetta;
    delete []barRreturn;


    calcF_D_FTransp(arrMtrx_dFGr_po_dX_Inv,arrTemp1,mDimX, arr_KsinsQ_per_1_rad);
    delete []arrTemp0;
    delete []arrTemp1;
    delete []arrTemp2;

     return true ;
}

//-----------------------------------------
// Вычисление кореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
bool  QLblSolver::calc_Ksins_Psi_Tetta(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kgps_Psi_Tetta_per_1_rad)
{

    memset(arr_Kgps_Psi_Tetta_per_1_rad, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();

    //double *parrArr_fi = new double [QntMeas * mDimX];
    double *parrArr_df_po_dPsiTetta = new double [QntMeas * mDimX*2];

    bool *barRreturn = new bool[QntMeas];

    memset(barRreturn, false,QntMeas * sizeof(bool));

    double arrPrirash[2] = {0.1,0.1};

                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrfi_and_arrdfi_po_dPsiTetta(arrX0
                                            ,i, arrPrirash,&parrArr_df_po_dPsiTetta[i * mDimX* 2])
               )
        {

            barRreturn[i] =  true;
        }

    }


 } // OMP (начало блока, который выполняется в нескольких потоках !

    for (int i = 0; i < QntMeas; ++i)
    {
         if (barRreturn [i])
         {

            // delete []parrArr_fi;
             delete []parrArr_df_po_dPsiTetta;
             delete []barRreturn;

        return false;
        }
    }


    double *arrTemp0 = new double [mDimX *mDimX];
    double *arrTemp1 = new double [mDimX *mDimX];
    memset(arrTemp1, 0, sizeof(double) *mDimX *mDimX);


    for (int i = 0; i < QntMeas; ++i)
    {
        calcF_Mult_FTransp(&parrArr_df_po_dPsiTetta[i * mDimX* 2],mDimX, 2, arrTemp0);
        MtrxSumMatrx(arrTemp1, arrTemp0,mDimX, mDimX, arrTemp1) ;

    }

   //delete []parrArr_fi;
    delete []parrArr_df_po_dPsiTetta;
    delete []barRreturn;


    calcF_D_FTransp(arrMtrx_dFGr_po_dX_Inv,arrTemp1,mDimX, arr_Kgps_Psi_Tetta_per_1_rad);
    delete []arrTemp0;
    delete []arrTemp1;


     return true ;
}

//--------------------------
bool  QLblSolver::calc_dX_po_dProfile(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dProfile)
{
    memset(arr_dX_po_dProfile, 0, mDimX * sizeof(double));
    double * parrMtrx_dFGr_po_dProfile = new double [mDimX];
    memset (parrMtrx_dFGr_po_dProfile, 0, mDimX* sizeof(double));

    calc_dFGr_po_dProfile(arrX0,  parrMtrx_dFGr_po_dProfile);
    //

    MtrxMultMatrx( arrMtrx_dFGr_po_dX_Inv,mDimX, mDimX, parrMtrx_dFGr_po_dProfile,1, arr_dX_po_dProfile) ;
    MatrxMultScalar(arr_dX_po_dProfile, 1, mDimX, -1.,arr_dX_po_dProfile);

    delete []parrMtrx_dFGr_po_dProfile;


     return true ;
}
//--------------------------------------------

bool  QLblSolver::calc_arrfi_and_arrdfi_po_dProfile(double *arrX
                       ,const int NUmBigMeasure,double *arrfi,double *arrdfi_po_dProfile)
{
   double arrNeviaz[8]= {0.};
    // !ОТЛАДКА
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        return false;
    }
     // вектор приращений
    double valPrirashen= 1.;


    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX ];
    double *arrNeviazCur = new double [mDimY];

        QLblSolver LblSolverTemp = *this;

       LblSolverTemp.mtblEstPrfl.addValue(valPrirashen);
      if (!LblSolverTemp.calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {

          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dProfile);

 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;

    return true;
}
//----------------------------------------
bool  QLblSolver::calc_KProfile(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KProfile_per_1_m)
{

    memset(arr_KProfile_per_1_m, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();
    int numProfileParts= mtblEstPrfl.mNumCols;
    double *parr_KPart = new double [numProfileParts * mDimX* mDimX];

    bool breturn = true;
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)

    for (int i = 0; i < numProfileParts; ++i)
    {
        if (!calc_K_ProfilePart(arrX0,i, arrMtrx_dFGr_po_dX_Inv,&parr_KPart[i *mDimX* mDimX])
               )
        {
            breturn =  false;
        }
    }

 } // OMP (начало блока, который выполняется в нескольких потоках !

    if(!breturn)
    {
        delete []parr_KPart;

        return false;
    }

    for (int i = 0; i < numProfileParts; ++i)
    {
        MtrxSumMatrx(arr_KProfile_per_1_m, &parr_KPart[i *mDimX* mDimX],mDimX, mDimX, arr_KProfile_per_1_m) ;
    }

     delete []parr_KPart;

     return true ;
}

//-------------------------------

bool  QLblSolver::calc_K_ProfilePart(double *arrX0, int numProfilePart,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KPart)
{
    memset(arr_KPart, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();

    double *parrArr_df_po_dPart = new double [mDimX];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));


    double *arrSum = new double [mDimX * mDimX];
    double *arrTemp = new double [mDimX * mDimX];
    memset (arrSum, 0, mDimX * mDimX * sizeof(double));
    for (int i = 0; i < QntMeas; ++i)
    {
        if (!calc_arrdfi_po_dProfilePart(arrX0,numProfilePart
                                            ,i,parrArr_df_po_dPart)
               )
        {

            barRreturn[i] =  true;
            break;
        }
        calcF_Mult_FTransp(parrArr_df_po_dPart,mDimX, 1, arrTemp);
        MtrxSumMatrx(arrSum, arrTemp,mDimX, mDimX, arrSum) ;

    }
    for (int i = 0; i < QntMeas; ++i)
    {
         if (barRreturn [i])
         {

             delete []parrArr_df_po_dPart;
             delete []arrSum;
             delete []arrTemp;
             delete []barRreturn;

        return false;
        }
    }

    calcF_D_FTransp(arrMtrx_dFGr_po_dX_Inv,arrSum,mDimX, arr_KPart);

    delete []parrArr_df_po_dPart;
    delete []arrSum;
    delete []arrTemp;
    delete []barRreturn;

     return true ;
}

//-------------------------------------------------
bool  QLblSolver::calc_arrdfi_po_dProfilePart(double *arrX, const int numProfilePart
                       ,const int NUmBigMeasure,double *arrdfi_po_dProfilePart)
{
   double arrNeviaz[8]= {0.};
   double *arrfi = new double [mDimX];
    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        return false;
    }
     // вектор приращений
    double valPrirashen = 1.;


    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [mDimX * mDimX];
    double *arrNeviazCur = new double [mDimY];


        QLblSolver LblSolverTemp = *this;

        LblSolverTemp.mtblEstPrfl.mparrVal[numProfilePart] += valPrirashen;


      if (!LblSolverTemp.calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {

          delete []arrfiCur;
          delete []arrfi;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dProfilePart);

 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrfi;

    return true;

}
//-----------------------------
void QLblSolver::calcArrPenalty(const double *arrX0,const double valPenalt, double *arrPenalt)
{
  memset(arrPenalt, 0, sizeof(double) *mDimX) ;
  double del1 = arrX0[2] + 15.;
  if (del1 < 0.)
  {
    arrPenalt[2]+=   valPenalt * del1 * del1;
  }
  double del2 = arrX0[2] + 3.;
  if (del2 > 2.)
  {
    arrPenalt[2]+=   valPenalt * del2 * del2;
  }
}
