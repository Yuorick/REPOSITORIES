#include "mythread.h"
#include <QDebug>
#include "LblSolver.h"
#include "Solver2D_Angs.h"
#include "Solver3D_Angs.h"
#include "BigMeasure.h"
#include "MatrixProccess.h"
#include<math.h>


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
                                ,mHBeaconEst, 0.015);

        double arrX0[8] ={0.},arrXRez [8] = {0.}, valNeviaz0 = -1.;
        for (int i = 0; i < 8; ++i)
        {
           arrX0[i] =  mvctX0.at(i);
        }
        QPosSolver *pPosSolver = &LblSolver;
        double sigmaTreshold = 3.;

        QVector<int> vctGoodZamersNum;
        pPosSolver->repaireZamerArray_(arrX0, sigmaTreshold,0, &vctGoodZamersNum);


        if (vctGoodZamersNum.size() <= 10)
        {
            delete []parrMeas;
            emit calcfinished(MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_LBL);
        }
        // построение облаков точек
           double arrMean[3] ={0.},arrDisp[3] ={0.};
           QVector<int> vctNumGood;
           //

           QString qstrOutPutFold0;
           wchar_t wchFolg0[400] = {0};
           if (mqstrOutPutFold.length() >3)
           {
            qstrOutPutFold0 = mqstrOutPutFold + QString("//CloudAfterRepaire_LBL").arg(0);
            qstrOutPutFold0.toWCharArray(wchFolg0);

            wchFolg0[qstrOutPutFold0.length()]=0;
            _wmkdir(wchFolg0);
           }
           //


            double arrBeaconPos[3] = {0.},arrAntPosParams [6] = {0.};
            arrBeaconPos[0] = mvctX0.at(3);
            arrBeaconPos[1] = mvctX0.at(4);
            arrBeaconPos[2] = mHBeaconEst;
            memcpy(arrAntPosParams,arrX0, 3 * sizeof(double) );
            memcpy(&arrAntPosParams[3],&arrX0[5], 3 * sizeof(double) );

           QPntCloud :: createPictFiles(wchFolg0,pPosSolver->mVectBigMeasures,mtblEstPrfl
                                       ,arrBeaconPos,arrAntPosParams, mTypeOfSolverTask
                                        ,arrMean, arrDisp,&vctNumGood,mType_of_Output_File
                                      );
           QVector <double> vctMean(3);
           QVector <double> vctDisp(3);
           for (int i =0;i< 3; ++i)
           {
             vctMean.replace(i,arrMean[i]) ;
             vctDisp.replace(i,arrDisp[i]) ;
           }
           kleaningEnded(1, vctGoodZamersNum,vctMean, vctDisp);
       // !0



        mvctMeas = pPosSolver->mVectBigMeasures;

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

      // выкинуть в mainwindow
      // построение облаков точек после решения LBL
      memset(wchFolg0,0, 400 * sizeof(wchar_t));

      if (mqstrOutPutFold.length() >3)
      {
       qstrOutPutFold0 = mqstrOutPutFold + QString("//CloudAfterSolv_LBL");
       qstrOutPutFold0.toWCharArray(wchFolg0);

       wchFolg0[qstrOutPutFold0.length()]=0;
       _wmkdir(wchFolg0);
      }
      /////

          arrBeaconPos[0] = arrXRez[3];
          arrBeaconPos[1] = arrXRez[4];
          arrBeaconPos[2] = mHBeaconEst;
          memcpy(arrAntPosParams,arrXRez, 3 * sizeof(double) );
          memcpy(&arrAntPosParams[3],&arrX0[5], 3 * sizeof(double) );

         QPntCloud :: createPictFiles(wchFolg0,pPosSolver->mVectBigMeasures,mtblEstPrfl
                                     ,arrBeaconPos,arrAntPosParams, mTypeOfSolverTask
                                      ,arrMean, arrDisp,&vctNumGood,mType_of_Output_File
                                    );

         for (int i =0;i< 3; ++i)
         {
           vctMean.replace(i,arrMean[i]) ;
           vctDisp.replace(i,arrDisp[i]) ;
         }
         kleaningEnded(2, vctGoodZamersNum,vctMean, vctDisp);
     // !



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
                 ,mHBeaconEst,0.002,arrXRez,arrSBeacon_GSK);

            pSolver2D_Angs = &Solver2D_Angs;
            pPosSolver = &Solver2D_Angs;

            break;

        case USBL_3D:
        Solver3D_Angs = QSolver3D_Angs(parrMeas, mvctMeas.size(),mtblEstPrfl
             ,mHBeaconEst,0.002,arrXRez,arrSBeacon_GSK);

        pSolver2D_Angs = &Solver3D_Angs;
        pPosSolver = &Solver3D_Angs;

            break;

        default:

            break;
        }


    // чистка массива замеров по курс углу
    double valTreshold = 15./ 180. * M_PI;
   pPosSolver->fncRoughCleaningMeasures(&arrX0[5],valTreshold, 0, &vctGoodZamersNum);
    sigmaTreshold = 3.;
    pPosSolver->repaireZamerArray_(&arrX0[5], sigmaTreshold,0, &vctGoodZamersNum);





   if (mTypeOfSolverTask == USBL_3D)
   {
    // чистка массива замеров по углу места
       pPosSolver->fncRoughCleaningMeasures(&arrX0[5],valTreshold, 1,&vctGoodZamersNum);
       pPosSolver->repaireZamerArray_(&arrX0[5], sigmaTreshold,1, &vctGoodZamersNum);

    //
     }

   if (vctGoodZamersNum.size() <= 10)
   {
       delete []parrMeas;
       emit calcfinished(MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_USBL);
   }

   // выкинуть в mainwindow
   // построение облаков точек после после чистки по невязке перед USBL



       memset(wchFolg0,0, 400 * sizeof(wchar_t));

       if (mqstrOutPutFold.length() >3)
       {
        qstrOutPutFold0 = mqstrOutPutFold + QString("//CloudAfterCleaning_USBL");
        qstrOutPutFold0.toWCharArray(wchFolg0);

        wchFolg0[qstrOutPutFold0.length()]=0;
        _wmkdir(wchFolg0);
       }
       arrBeaconPos[0] = arrXRez[3];
       arrBeaconPos[1] = arrXRez[4];
       arrBeaconPos[2] = mHBeaconEst;
       memcpy(arrAntPosParams,arrXRez, 3 * sizeof(double) );
       memcpy(&arrAntPosParams[3],&arrX0[5], 3 * sizeof(double) );

      QPntCloud :: createPictFiles(wchFolg0,pPosSolver->mVectBigMeasures,mtblEstPrfl
                                  ,arrBeaconPos,arrAntPosParams, mTypeOfSolverTask
                                   ,arrMean, arrDisp,&vctNumGood,mType_of_Output_File
                                 );

      for (int i =0;i< 3; ++i)
      {
        vctMean.replace(i,arrMean[i]) ;
        vctDisp.replace(i,arrDisp[i]) ;
      }
      kleaningEnded(3, vctGoodZamersNum,vctMean, vctDisp);
  // !


   mvctMeas = pPosSolver->mVectBigMeasures;
       double arrMtrx_dFGr_po_dX_Inv1[100] = {0.};
        if (mALGOR_TYPE == PEREBOR)
        {
            THREAD_RESULT = estimateAngs(pSolver2D_Angs,  &arrX0[5]
            ,numPart, mMaxQuantIter_USBL,&arrXRez[5]
            ,valNeviaz0, arrMtrx_dFGr_po_dX_Inv);
        }
        if (mALGOR_TYPE == NUTON)
        {
            THREAD_RESULT = estimateParams(pPosSolver,  &arrX0[5]
            ,numPart ,mMaxQuantIter_USBL, mUSBLcoeff, &arrXRez[5]
            , valNeviaz0, arrMtrx_dFGr_po_dX_Inv1);
        }

        // выкинуть в mainwindow
        // построение облаков точек после после после решения USBL


            memset(wchFolg0,0, 400 * sizeof(wchar_t));

            if (mqstrOutPutFold.length() >3)
            {
             qstrOutPutFold0 = mqstrOutPutFold + QString("//CloudAfterSolving_USBL");
             qstrOutPutFold0.toWCharArray(wchFolg0);

             wchFolg0[qstrOutPutFold0.length()]=0;
             _wmkdir(wchFolg0);
            }

            arrBeaconPos[0] = arrXRez[3];
            arrBeaconPos[1] = arrXRez[4];
            arrBeaconPos[2] = mHBeaconEst;
            memcpy(arrAntPosParams,arrXRez, 3 * sizeof(double) );
            memcpy(&arrAntPosParams[3],&arrXRez[5], 3 * sizeof(double) );

           QPntCloud :: createPictFiles(wchFolg0,pPosSolver->mVectBigMeasures,mtblEstPrfl
                                       ,arrBeaconPos,arrAntPosParams, mTypeOfSolverTask
                                        ,arrMean, arrDisp,&vctNumGood,mType_of_Output_File
                                      );

           for (int i =0;i< 3; ++i)
           {
             vctMean.replace(i,arrMean[i]) ;
             vctDisp.replace(i,arrDisp[i]) ;
           }
           kleaningEnded(4, vctGoodZamersNum,vctMean, vctDisp);
       // !

    delete []parrMeas;


    // отправляем в интерфейс результат выполнения задачи
    emit calcfinished(THREAD_RESULT);
}


//--------------------------------
void MyThread::setParams(QVector <QBigMeasure> vctMeas, TYPE_OF_SOLVER_TASK TypeOfSolverTask
                        ,QVector <double> vctX0,TTable_1D tblEstPrfl, double *arrGPSAprioriPosParams
                       ,double HBeaconEst,int quantIter_LBL, int quantIter_USBL
                         ,double Gape,double AngStep, ALGOR_TYPE nALGOR_TYPE,double LBLcoeff,double USBLcoeff
                         , TYPE_OF_OUTPUT_FILE Type_of_Output_File, QString qstrOutPutFold)
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
   mType_of_Output_File = Type_of_Output_File;
   mqstrOutPutFold = qstrOutPutFold;
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
   double valGape_bet = mGape/180. * M_PI;
   double valGape_eps = ((mGape > 45.)?45.:mGape)/180. * M_PI;
   double valGape_alf = 16./180.* M_PI;

   double in_v_storonu_bet = valGape_bet/arg_step;
   int in_bet = 2* in_v_storonu_bet +1;

   double in_v_storonu_eps = valGape_eps/arg_step;
   int in_eps = 2* in_v_storonu_eps +1;

   double in_v_storonu_alf = valGape_alf/arg_step;
   int in_alf = 2* in_v_storonu_alf +1;


   double bet_min = arrX0[0] - ((double)in_v_storonu_bet) * arg_step;
   double eps_min = arrX0[1] - ((double)in_v_storonu_eps) * arg_step;
   double alf_min = arrX0[2] - ((double)in_v_storonu_alf) * arg_step;



   double valNevCur = 0.;

   int iChetchik = 0;
   for (int i =0; i < in_bet; ++i)
       for (int j =0; j < in_eps; ++j)
           for (int k =0; k < in_alf; ++k)
   {
      double arrXcur[3] = {0.};
      arrXcur[0] = bet_min + ((double)i) * arg_step;
      arrXcur[1] = eps_min + ((double)j) * arg_step;
      arrXcur[2] = alf_min + ((double)k) * arg_step;

      valNevCur = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXcur, quantGoodZamers0
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

   int nSt = in_bet * in_eps * in_alf;
emit progress(nSt, vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);
return MYTHREAD_RESULT_OK;

}



