#include "mythread.h"
#include <QDebug>
#include "LblSolver.h"
#include "Solver2D_Angs.h"
#include "Solver3D_Angs.h"
#include "BigMeasure.h"
#include "MatrixProccess.h"
#include "CoordSystTrsf.h"
#include <math.h>
#include <float.h>

MyThread::MyThread(QObject *parent) :
    QThread(parent),
  userBreakRequested(false)
 // ,maxStep(1),
 // errorThreshold(0.1)
{    
}

void MyThread::run()
{
    qDebug() << "MyThread::run(), in";

    int result = MYTHREAD_RESULT_MAX_STEP; // результат по умолчанию
    int icur = 0;

    int  numPart = 1;



        // !следующую строку следует заменить на реальные вычисления!
        //QThread::msleep(100); // этот вызов здесь нужен для имитации длительности вычислений
      //****************************
        //const int QuantMeas = ;
         QBigMeasure *parrMeas = (QBigMeasure *)malloc((mvctMeas.size()) * sizeof(QBigMeasure) );
         for (int i = 0; i < mvctMeas.size(); ++i )
         {
            parrMeas[i] =  mvctMeas.at(i);
         }
      QLblSolver LblSolver(parrMeas, mvctMeas.size()    ,mtblEstPrfl
                                ,marrGPSAprioriPosParams,mHBeaconEst, 0.015);



        double arrX0[8] ={0.},arrXRez [8] = {0.}, valNeviaz0 = -1.;
        for (int i = 0; i < 8; ++i)
        {
           arrX0[i] =  mvctX0.at(i);
        }
        QPosSolver *pPosSolver = &LblSolver;
        double sigmaTreshold = 2.5;
        pPosSolver->repaireZamerArray_(arrX0, sigmaTreshold);
        mvctMeas = pPosSolver->mVectBigMeasures;

       // mvctMeas = pPosSolver->mVectBigMeasures;


        double arrMtrx_dFGr_po_dX_Inv[100] = {0.};

        MYTHREAD_RESULT THREAD_RESULT = estimateParams(pPosSolver,  arrX0
             ,1,mMaxQuantIter_LBL, mLBLcoeff,arrXRez, valNeviaz0, arrMtrx_dFGr_po_dX_Inv);

      switch (THREAD_RESULT)
      {

      case MYTHREAD_RESULT_MAX_STEP:   // вычисления провалены (достигли максимально допустимого кол-ва итераций)
      case MYTHREAD_RESULT_USER_CANCEL: // вычисления провалены (прервано пользователем)
      case MYTHREAD_RESULT_ERROR:
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
          break;

      default:
          break;
      }

      if (mTypeOfSolverTask == LBL)
      {
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
      }


//***************************************
      numPart = 2;

    qDebug() << "MyThread::run(), out";

    //*********************************
    double arrSBeacon_GSK[3] = {0.};
    arrSBeacon_GSK[0] = arrXRez[3];
    arrSBeacon_GSK[1] = arrXRez[4];
    arrSBeacon_GSK[2] = mHBeaconEst;

    QSolver3D_Angs Solver3D_Angs;
    QSolver2D_Angs Solver2D_Angs;
    QSolver2D_Angs *pSolver2D_Angs;


    switch(mTypeOfSolverTask )
        {

        case USBL_2D:
            Solver2D_Angs = QSolver2D_Angs(parrMeas, mvctMeas.size(),mtblEstPrfl
                 ,marrGPSAprioriPosParams,mHBeaconEst,0.009,arrXRez,arrSBeacon_GSK);

            pSolver2D_Angs = &Solver2D_Angs;
            pPosSolver = &Solver2D_Angs;

            break;

        case 2:
        Solver3D_Angs = QSolver3D_Angs(parrMeas, mvctMeas.size(),mtblEstPrfl
             ,marrGPSAprioriPosParams,mHBeaconEst,0.009,arrXRez,arrSBeacon_GSK);

        pSolver2D_Angs = &Solver3D_Angs;
        pPosSolver = &Solver3D_Angs;

            break;

        default:

            break;
        }
    //sigmaTreshold = 2.2;
   // pPosSolver->repaireZamerArray_(arrX0, sigmaTreshold);
    double valTreshold = 5./ 180. * M_PI;
    pPosSolver->fncRoughCleaningMeasures(arrX0,valTreshold);
    mvctMeas = pPosSolver->mVectBigMeasures;

       double arrMtrx_dFGr_po_dX_Inv1[100] = {0.};
        if (mALGOR_TYPE == PEREBOR)
        {
            THREAD_RESULT = estimateAngs(pSolver2D_Angs,  &arrX0[5]
            , numPart, mMaxQuantIter_USBL,&arrXRez[5]
            , valNeviaz0, arrMtrx_dFGr_po_dX_Inv);
        }
        if (mALGOR_TYPE == NUTON)
        {
            THREAD_RESULT = estimateParams(pPosSolver,  &arrX0[5]
            ,numPart ,mMaxQuantIter_USBL, mUSBLcoeff, &arrXRez[5]
            , valNeviaz0, arrMtrx_dFGr_po_dX_Inv1);
        }

    delete []parrMeas;


    // отправляем в интерфейс результат выполнения задачи
    emit calcfinished(THREAD_RESULT);
}


//--------------------------------
void MyThread::setParams(QVector <QBigMeasure> vctMeas, TYPE_OF_SOLVER_TASK TypeOfSolverTask
                        ,QVector <double> vctX0,TTable_1D tblEstPrfl, double *arrGPSAprioriPosParams
                       ,double HBeaconEst,int quantIter_LBL, int quantIter_USBL
                         ,double Gape,double AngStep, ALGOR_TYPE nALGOR_TYPE,double LBLcoeff,double USBLcoeff)
{
   mvctMeas = vctMeas;
   mTypeOfSolverTask = TypeOfSolverTask;
   mvctX0 = vctX0;
   mtblEstPrfl = tblEstPrfl;
   memcpy(marrGPSAprioriPosParams, arrGPSAprioriPosParams, 3 * sizeof(double));
   mHBeaconEst = HBeaconEst;

   mMaxQuantIter_LBL = quantIter_LBL;
   mMaxQuantIter_USBL = quantIter_USBL;
   mGape =Gape;
   mAngStep = AngStep;
   mALGOR_TYPE = nALGOR_TYPE;
   mLBLcoeff = LBLcoeff;
   mUSBLcoeff = USBLcoeff;
}
//--------------------------------

// вычисление неизвестного вектора позиционирования методом наименьших квадратов
//arrX[0],arrX[1],arrX[2] - вектор параллакса антенны ПСК
//arrX[3],arrX[4],arrX[5] - вектор координат маяка в ГСК
//arrX[6],arrX[7],arrX[8] - вектор углов ориентации антенны в ПСК
// INPUT:
//arrX0[mDimX] - начальное приближение
// OUTPUT:
//arrXEst[mDimX] - решение
// valNeviaz0 - невязка
MYTHREAD_RESULT MyThread::estimateParams(QPosSolver *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter, double coeff,double *arrXRez
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{

    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));

   double *arrMean = new double [mDimY];
   double *arrDisp = new double [mDimY];
   double *arrNeviazSquare = new double [mDimY];


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                , arrMean, arrDisp, arrNeviazSquare);
   QVector<double> vctMean(mDimY);
   QVector<double> vctDisp(mDimY);
   QVector<double> vctNeviazSquare(mDimY);
   for (int i =0; i < mDimY;++i)
   {
     vctMean.replace(i,arrMean[i]) ;
     vctDisp.replace(i,arrDisp[i]) ;
     vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
   }
   ///

   // 3. итреац процесс
  double *arrXCur = new double[mDimX];
  memcpy(arrXCur, arrXRez, mDimX * sizeof(double));

  double valNeviaz1=0.;
  double *arrDel = new double[mDimX];

  memset(arrDel, 0, mDimX* sizeof(double));

  bool breturn = false;

   for (int i =0; i < MaxQuantIter; ++i)
   {
       if(!(pPosSolver->doOneIteration(arrXRez, arrDel,arrMtrx_dFGr_po_dX_Inv)))
       {
       delete []arrXCur;
       delete []arrDel;
       return MYTHREAD_RESULT_ERROR;

       }

       double del1 = NormVect(arrDel, mDimX);

       MatrxMultScalar( arrDel, mDimX, 1, coeff, arrDel);

       MtrxMinusMatrx(arrXRez, arrDel,mDimX, 1, arrXCur);

       memcpy(arrXRez, arrXCur, mDimX * sizeof(double));


      pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                   , arrMean, arrDisp, arrNeviazSquare);

      for (int i = 0; i < mDimX; ++i)
      {
          mvctX0.replace((numPart -1) *5 +i, arrXRez[i]);

      }
      for (int i =0; i < mDimY;++i)
      {
        vctMean.replace(i,arrMean[i]) ;
        vctDisp.replace(i,arrDisp[i]) ;
        vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
      }
       if (del1 < pPosSolver->mToler)
       {
           delete []arrXCur;
           delete []arrDel;
           delete[]arrMean ;
           delete[]arrDisp ;
           delete[]arrNeviazSquare ;

        emit progress(i , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);

       return MYTHREAD_RESULT_OK;

       }

      if(userBreakRequested)
      {
          delete []arrXCur;
          delete []arrDel;
          delete[]arrMean ;
          delete[]arrDisp ;
          delete[]arrNeviazSquare ;
          return MYTHREAD_RESULT_USER_CANCEL;
      }

       emit progress(i , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);
   }



   delete []arrXCur;
   delete []arrDel;
   delete[]arrMean ;
   delete[]arrDisp ;
   delete[]arrNeviazSquare ;

    return MYTHREAD_RESULT_MAX_STEP;
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
MYTHREAD_RESULT MyThread::estimateAngs(QSolver2D_Angs *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter,double *arrXRez
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{


    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;

    double *arrMean = new double [mDimY];
    double *arrDisp = new double [mDimY];
    double *arrNeviazSquare = new double [mDimY];
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                        , arrMean, arrDisp, arrNeviazSquare);
   QVector<double> vctMean(mDimY);
   QVector<double> vctDisp(mDimY);
   QVector<double> vctNeviazSquare(mDimY);
   for (int i =0; i < mDimY;++i)
   {
     vctMean.replace(i,arrMean[i]) ;
     vctDisp.replace(i,arrDisp[i]) ;
     vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
   }
   ///
   // перебор 3-ной цикл
   double arg_step = mAngStep/ 180. * M_PI;
   double valGape = mGape/180. * M_PI;
   double in_v_storonu = valGape/arg_step;
   int in = 2* in_v_storonu +1;
   double bet_min = arrX0[0] - ((double)in_v_storonu) * arg_step;
   double eps_min = arrX0[1] - ((double)in_v_storonu) * arg_step;
   double alf_min = arrX0[2] - ((double)in_v_storonu) * arg_step;

   int nSt = in * in * in;

   double valNevCur = 0.;

   int iChetchik = 0;
   for (int i =0; i < in; ++i)
       for (int j =0; j < in; ++j)
           for (int k =0; k < in; ++k)
   {
      double arrXcur[3] = {0.};
      arrXcur[0] = bet_min + ((double)i) * arg_step;
      arrXcur[1] = eps_min + ((double)j) * arg_step;
      arrXcur[2] = alf_min + ((double)k) * arg_step;

      valNevCur = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                          , arrMean, arrDisp, arrNeviazSquare);
      if (valNevCur < valNeviaz0)
      {
         valNeviaz0 =  valNevCur;
         memcpy(arrXRez, arrXcur,3 * sizeof(double));
         for (int i = 0; i < mDimX; ++i)
         {
             mvctX0.replace((numPart -1) *5 +i, arrXRez[i]);

         }
         for (int i =0; i < mDimY;++i)
         {
           vctMean.replace(i,arrMean[i]) ;
           vctDisp.replace(i,arrDisp[i]) ;
           vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
         }

      }
      ++iChetchik;
      if ((iChetchik % 50) == 0)
      {
         emit progress(iChetchik , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);

          if(userBreakRequested)
          {
              delete[]arrMean ;
              delete[]arrDisp ;
              delete[]arrNeviazSquare ;

              return MYTHREAD_RESULT_USER_CANCEL;
          }
      }

   }


   delete[]arrMean ;
   delete[]arrDisp ;
   delete[]arrNeviazSquare ;

emit progress(nSt, vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);
return MYTHREAD_RESULT_OK;

}

/*
void MyThread::run()
{
    qDebug() << "MyThread::run(), in";

    int result = MYTHREAD_RESULT_MAX_STEP; // результат по умолчанию
    int icur = 0;

    int  numPart = 1;



        // !следующую строку следует заменить на реальные вычисления!
        //QThread::msleep(100); // этот вызов здесь нужен для имитации длительности вычислений
      //****************************
        const int QuantMeas = mvctMeas.size();
         QBigMeasure *parrMeas = (QBigMeasure *)malloc(QuantMeas * sizeof(QBigMeasure) );
         for (int i = 0; i <QuantMeas; ++i )
         {
            parrMeas[i] =  mvctMeas.at(i);
         }
      QLblSolver LblSolver(parrMeas, QuantMeas    ,mtblEstPrfl
                                ,marrGPSAprioriPosParams,mHBeaconEst, 0.015);

        double arrX0[8] ={0.},arrXRez [8] = {0.}, valNeviaz0 = -1.;
        for (int i = 0; i < 8; ++i)
        {
           arrX0[i] =  mvctX0.at(i);
        }

        QPosSolver *pPosSolver = &LblSolver;
        double arrMtrx_dFGr_po_dX_Inv[100] = {0.};

        MYTHREAD_RESULT THREAD_RESULT = estimateParams(pPosSolver,  arrX0
             ,1,mMaxQuantIter_LBL, arrXRez, 0.25, valNeviaz0, arrMtrx_dFGr_po_dX_Inv);



      switch (THREAD_RESULT)
      {

      case MYTHREAD_RESULT_MAX_STEP:   // вычисления провалены (достигли максимально допустимого кол-ва итераций)
      case MYTHREAD_RESULT_USER_CANCEL: // вычисления провалены (прервано пользователем)
      case MYTHREAD_RESULT_ERROR:
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
          break;

      default:
          break;
      }

      if (mTypeOfSolverTask == LBL)
      {
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
      }


//***************************************
      numPart = 2;

    qDebug() << "MyThread::run(), out";

    //*********************************
    double arrSBeacon_GSK[3] = {0.};
    arrSBeacon_GSK[0] = arrXRez[3];
    arrSBeacon_GSK[1] = arrXRez[4];
    arrSBeacon_GSK[2] = mHBeaconEst;


    QSolver3D_Angs Solver3D_Angs;
    QSolver2D_Angs Solver2D_Angs;
    QSolver2D_Angs *pSolver2D_Angs;


    switch(mTypeOfSolverTask )
        {

        case USBL_2D:
            Solver2D_Angs = QSolver2D_Angs(parrMeas, QuantMeas    ,mtblEstPrfl
                 ,marrGPSAprioriPosParams,mHBeaconEst,0.009,arrXRez,arrSBeacon_GSK);

            pSolver2D_Angs = &Solver2D_Angs;
            pPosSolver = &Solver2D_Angs;

            break;

        case 2:
        Solver3D_Angs = QSolver3D_Angs(parrMeas, QuantMeas    ,mtblEstPrfl
             ,marrGPSAprioriPosParams,mHBeaconEst,0.009,arrXRez,arrSBeacon_GSK);

        pSolver2D_Angs = &Solver3D_Angs;
        pPosSolver = &Solver3D_Angs;

            break;

        default:

            break;
        }


       double arrMtrx_dFGr_po_dX_Inv1[100] = {0.};
        if (mALGOR_TYPE == PEREBOR)
        {
            THREAD_RESULT = estimateAngs(pSolver2D_Angs,  &arrX0[5]
            , numPart, mMaxQuantIter_USBL,&arrXRez[5]
            , valNeviaz0, arrMtrx_dFGr_po_dX_Inv);
        }
        if (mALGOR_TYPE == NUTON)
        {
            THREAD_RESULT = estimateParams(pPosSolver,  &arrX0[5]
            ,numPart ,mMaxQuantIter_USBL, &arrXRez[5], mcoeff
            , valNeviaz0, arrMtrx_dFGr_po_dX_Inv1);
        }

    delete []parrMeas;

       // ОТЛАДКА
      // double matrPereh_PSK_V_KGSK_REZ[9] = {0.}, matrPereh_PSK_V_KGSK_0[9] = {0.};
      // QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&arrXRez[5],matrPereh_PSK_V_KGSK_REZ);
      // QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&arrX0[5],matrPereh_PSK_V_KGSK_0);
       // !ОТЛАДКА

    // отправляем в интерфейс результат выполнения задачи
    emit calcfinished(THREAD_RESULT);
}
//--------------------------------
void MyThread::setParams(QVector <QBigMeasure> vctMeas, TYPE_OF_SOLVER_TASK TypeOfSolverTask
                        ,QVector <double> vctX0,TTable_1D tblEstPrfl, double *arrGPSAprioriPosParams
                       ,double HBeaconEst,int quantIter_LBL, int quantIter_USBL
                         ,double Gape,double AngStep, ALGOR_TYPE nALGOR_TYPE,double coeff)


{
   mvctMeas = vctMeas;
   mTypeOfSolverTask = TypeOfSolverTask;
   mvctX0 = vctX0;
   mtblEstPrfl = tblEstPrfl;
   memcpy(marrGPSAprioriPosParams, arrGPSAprioriPosParams, 3 * sizeof(double));
   mHBeaconEst = HBeaconEst;

   mMaxQuantIter_LBL = quantIter_LBL;
   mMaxQuantIter_USBL = quantIter_USBL;
   mGape =Gape;
   mAngStep = AngStep;
   mALGOR_TYPE = nALGOR_TYPE;
   mcoeff = coeff;
}
//-----------------------------------------------------------------
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
MYTHREAD_RESULT MyThread::estimateParams(QPosSolver *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter,double *arrXRez, double coef
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{

    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcSumNeviazka(arrXRez, quantGoodZamers0);
   ///

   // 3. итреац процесс
  double *arrXCur = new double[mDimX];
  memcpy(arrXCur, arrXRez, mDimX * sizeof(double));

  double valNeviaz1=0.;
  double *arrDel = new double[mDimX];

  memset(arrDel, 0, mDimX* sizeof(double));

  bool breturn = false;

   for (int i =0; i < MaxQuantIter; ++i)
   {

       if(!(pPosSolver->doOneIteration(arrXRez, arrDel,arrMtrx_dFGr_po_dX_Inv)))
       {
       delete []arrXCur;
       delete []arrDel;
       return MYTHREAD_RESULT_ERROR;

       }

       int len = mDimX;
       if (mDimX == 8)
       {
       len = 5;
       }
       double del1 = NormVect(arrDel, len);
       double del2 = (mDimX == 8)?NormVect(&arrDel[5], 3):0.;


       MatrxMultScalar( arrDel, mDimX, 1, coef, arrDel);

       MtrxMinusMatrx(arrXRez, arrDel,mDimX, 1, arrXCur);

       memcpy(arrXRez, arrXCur, mDimX * sizeof(double));

       valNeviaz0 = pPosSolver->calcSumNeviazka(arrXRez, quantGoodZamers0);
       if ((del1 < pPosSolver->mToler)&&(del2 < 0.0005))
       {
           delete []arrXCur;
           delete []arrDel;

           for (int i = 0; i < mDimX; ++i)
           {
               mvctX0.replace((numPart -1) *5 +i, arrXRez[i]);

           }
        emit progress(i , sqrt(valNeviaz0/mvctMeas.size()/ (pPosSolver->mDimX)),numPart, true,mvctX0);

       return MYTHREAD_RESULT_OK;

       }

      if(userBreakRequested)
      {
          delete []arrXCur;
          delete []arrDel;
          return MYTHREAD_RESULT_USER_CANCEL;
      }

       emit progress(i , del1,numPart,false,mvctX0);
   }



   delete []arrXCur;
   delete []arrDel;

    return MYTHREAD_RESULT_MAX_STEP;
}
//------------------------------------
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
MYTHREAD_RESULT MyThread::estimateAngs(QSolver2D_Angs *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter,double *arrXRez
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{

    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcSumNeviazka(arrXRez, quantGoodZamers0);
   ///
   // перебор 3-ной цикл
   double arg_step = mAngStep/ 180. * M_PI;
   double valGape = mGape/180. * M_PI;
   double in_v_storonu = valGape/arg_step;
   int in = 2* in_v_storonu +1;
   double bet_min = arrX0[0] - ((double)in_v_storonu) * arg_step;
   double eps_min = arrX0[1] - ((double)in_v_storonu) * arg_step;
   double alf_min = arrX0[2] - ((double)in_v_storonu) * arg_step;

   int nSt = in * in * in;

   double valNevCur = 0.;

   int iChetchik = 0;
   for (int i =0; i < in; ++i)
       for (int j =0; j < in; ++j)
           for (int k =0; k < in; ++k)
   {
      double arrXcur[3] = {0.};
      arrXcur[0] = bet_min + ((double)i) * arg_step;
      arrXcur[1] = eps_min + ((double)j) * arg_step;
      arrXcur[2] = alf_min + ((double)k) * arg_step;

      valNevCur = pPosSolver->calcSumNeviazka(arrXcur, quantGoodZamers0);
      if (valNevCur < valNeviaz0)
      {
         valNeviaz0 =  valNevCur;
         memcpy(arrXRez, arrXcur,3 * sizeof(double));
         for (int i = 0; i < mDimX; ++i)
         {
             mvctX0.replace((numPart -1) *5 +i, arrXRez[i]);

         }

      }
      ++iChetchik;
      if ((iChetchik % 50) == 0)
      {
         emit progress(iChetchik , sqrt(valNeviaz0/mvctMeas.size()/3.),numPart, true,mvctX0);
          if(userBreakRequested)
          {

              return MYTHREAD_RESULT_USER_CANCEL;
          }
      }

   }


emit progress(nSt , sqrt(valNeviaz0/mvctMeas.size()),numPart, true,mvctX0);

return MYTHREAD_RESULT_OK;

}*/



